#=
function _unpack_local_solve!(m::GlobalOptimizer, opt::T) where T

    m._maximum_node_id = 0
    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())
    m._feasible_solution_found = m._result_status_code == MOI.FEASIBLE_POINT

    if MOI.get(opt, MOI.ResultCount()) > 0

        objective_value = MOI.get(opt, MOI.ObjectiveValue())

        # corrects for standard printing multiplier
        if _optimization_sense(m._input_problem) == MOI.MAX_SENSE
            objective_value *= -1.0
        end

        m._global_lower_bound = objective_value
        m._global_upper_bound = objective_value
        m._objective_value = objective_value
        m._best_upper_value = objective_value
        m._solution_value = objective_value
    end

    map!(x -> MOI.get(opt, MOI.VariablePrimal(), x), m._solution, m._relaxed_variable_index)

    m._run_time = time() - m._start_time

    return nothing
end

function add_sv_or_aff_obj!(m::GlobalOptimizer, opt::T) where T
    if _objective_type(_input_problem(m)) == SINGLE_VARIABLE
        MOI.set(opt, MOI.ObjectiveFunction{SV}(), m._input_problem._objective_sv)
    elseif _objective_type(_input_problem(m)) == SCALAR_AFFINE
        MOI.set(opt, MOI.ObjectiveFunction{SAF}(), m._input_problem._objective)
    end
    return nothing
end
=#

function add_variables(m::GlobalOptimizer, optimizer::T, variable_number::Int) where T

    variable_index = fill(VI(1), variable_number)
    for i = 1:variable_number
        @inbounds variable_index[i] = MOI.add_variable(optimizer)
        single_variable = SV(@inbounds variable_index[i])
        v_info = @inbounds m._working_problem._variable_info[i]
        if _is_fixed(v_info)
            MOI.add_constraint(optimizer, single_variable, _ET(vi))
        elseif _is_real_interval(vi)
            MOI.add_constraint(optimizer, single_variable, _IT(vi))
        elseif _is_greater_than(vi)
            MOI.add_constraint(optimizer, single_variable, _GT(vi))
        elseif _is_less_than(vi)
            MOI.add_constraint(optimizer, single_variable, _LT(vi))
        elseif _is_zero_one(vi)
            MOI.add_constraint(optimizer, single_variable, _ZO(vi))
        elseif _is_int_interval(vi)
            MOI.add_constraint(optimizer, single_variable, _Int(vi))
        end
    end

    return variable_index
end

"""
$(TYPEDSIGNATURES)

Loads variables, linear constraints, and empty storage for first nlp and
quadratic cut.
"""
function load_relaxed_problem!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    opt = m.relaxed_optimizer

    # add variables and indices and constraints
    wp = m._working_problem
    branch_num = 0

    variable_count = wp._variable_num
    for i = 1:variable_count

        vinfo = wp._variable_info[i]
        is_branch = m._branch_variables[i]
        is_branch && (branch_num += 1)

        var_indx = MOI.add_variable(opt)
        var = SV(var_indx)
        push!(m._relaxed_variable_index, var_indx)

        if vinfo.is_integer
        elseif vinfo.is_fixed
            ci_et = MOI.add_constraint(opt, var, ET(vinfo.lower_bound))
            is_branch && push!(m._relaxed_variable_eq, (ci_sv_et, branch_num))
        else
            if vinfo.has_lower_bound
                ci_gt = MOI.add_constraint(opt, var, GT(vinfo.lower_bound))
                is_branch && push!(m._relaxed_variable_gt, (ci_gt, branch_num))
            end
            if vinfo.has_upper_bound
                ci_lt = MOI.add_constraint(opt, var, LT(vinfo.upper_bound))
                is_branch && push!(m._relaxed_variable_lt, (ci_lt, branch_num))
            end
        end
    end

    # set node index to single variable constraint index maps
    m._branch_variable_num = branch_num
    m._node_to_sv_leq_ci = fill(CI{SV,LT}(-1), branch_num)
    m._node_to_sv_geq_ci = fill(CI{SV,GT}(-1), branch_num)
    for v in m._relaxed_variable_lt
        m._node_to_sv_leq_ci[v[2]] = v[1]
    end
    for v in m._relaxed_variable_gt
        m._node_to_sv_geq_ci[v[2]] = v[1]
    end

    # add linear constraints (ASSUME NOT BRIDGED... TODO: FIX LATER)
    _add_linear_constraints!(opt, m._input_problem)
    MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    return
end

function _mod_vinfo!(wp, ip, ci::CI{SV,S}) where S
    sv = MOI.get(ip, MOI.ConstraintFunction(), ci)
    sv_set = MOI.get(ip, MOI.ConstraintSet(), ci)
    wp._variable_info[sv.value] = VariableInfo(wp._variable_info[sv.value], sv_set)
    return nothing
end

function _add_decision_variables!(wp, ip::InputModel{T}) where {T<:Real}
    imodel = ip._input_model
    ip_variable_num = MOI.get(imodel, MOI.NumberOfVariables())
    wp._variable_num = ip_variable_num
    wp._variable_info = fill(VariableInfo{T}(), ip_variable_num)
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, LT{T}}()))
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, GT{T}}()))
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, ET{T}}()))
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, ZO}()))
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, MOI.Integer}()))
    any(isempty, wp._variable_info) && error("Variable bounds lead to infeasible problem.") # TODO: Message + Model Issue Return Code
    return nothing
end

function _add_constraint!(wp, ip::InputModel{T}) where {T<:Real}
    imodel = ip._input_model
    ip_variable_num = MOI.get(imodel, MOI.NumberOfVariables())
    wp._variable_num = ip_variable_num
    wp._variable_info = fill(VariableInfo{T}(), ip_variable_num)
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, LT{T}}()))
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, GT{T}}()))
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, ET{T}}()))
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, ZO}()))
    foreach(x -> _mod_vinfo!(wp, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SV, MOI.Integer}()))
    any(isempty, wp._variable_info) && error("Variable bounds lead to infeasible problem.") # TODO: Message + Model Issue Return Code
    return nothing
end

for Q in (LT,GT)
    @eval function _add_constraint!(wp::ParsedProblem, ip::InputModel{T}, ci::CI{SAF{T},$Q{T}}) where {T<:Real}
        imodel = ip._input_model
        f = MOI.get(imodel, MOI.ConstraintFunction(), ci)
        s = MOI.get(imodel, MOI.ConstraintSet(), ci)
        push!(wp._saf_leq, AffineFunctionIneq(f, s))
        return nothing
    end
    @eval function _add_constraint!(wp::ParsedProblem, ip::InputModel{T}, ci::CI{SQF{T},$Q{T}}) where {T<:Real}
        imodel = ip._input_model
        f = MOI.get(imodel, MOI.ConstraintFunction(), ci)
        s = MOI.get(imodel, MOI.ConstraintSet(), ci)
        push!(wp._sqf_leq, BufferedQuadraticIneq(f, s))
        return nothing
    end
end
function _add_constraint!(wp::ParsedProblem, ip::InputModel{T}, ci::CI{SAF{T},ET{T}}) where {T<:Real}
    imodel = ip._input_model
    f = MOI.get(imodel, MOI.ConstraintFunction(), ci)
    s = MOI.get(imodel, MOI.ConstraintSet(), ci)
    push!(wp._saf_eq, AffineFunctionEq(f, s))
    return nothing
end
function _add_constraint!(wp::ParsedProblem, ip::InputModel{T}, ci::CI{SQF{T},ET{T}}) where {T<:Real}
    imodel = ip._input_model
    f = MOI.get(imodel, MOI.ConstraintFunction(), ci)
    s = MOI.get(imodel, MOI.ConstraintSet(), ci)
    push!(wp._sqf_eq, BufferedQuadraticEq(f, s))
    return nothing
end

function _add_constraint!(wp::ParsedProblem, ip::InputModel{T}, fs::CI{VECVAR,SOC_CONE}) where {T<:Real}
    imodel = ip._input_model
    f = MOI.get(imodel, MOI.ConstraintFunction(), ci)
    s = MOI.get(imodel, MOI.ConstraintSet(), ci)
    first_variable_loc = f.variables[1].value
    prior_lbnd = wp._variable_info[first_variable_loc].lower_bound
    wp._variable_info[first_variable_loc].lower_bound = max(prior_lbnd, 0.0)
    push!(wp._conic_second_order, BufferedSOC(f, s))
    return nothing
end

function _add_linear_constraints!(m::S, ip::InputModel{T}) where {S,T<:Real}
    imodel = ip._input_model
    foreach(x -> _add_constraint!(m, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SAF{T}, LT{T}}()))
    foreach(x -> _add_constraint!(m, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SAF{T}, GT{T}}()))
    foreach(x -> _add_constraint!(m, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SAF{T}, ET{T}}()))
    return nothing
end

function _add_quadratic_constraints!(m::S, ip::InputModel{T}) where {S,T<:Real}
    imodel = ip._input_model
    foreach(x -> _add_constraint!(m, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SQF{T}, LT{T}}()))
    foreach(x -> _add_constraint!(m, imodel, x), MOI.get(imodel, MOI.ListOfConstraintIndices{SQF{T}, GT{T}}()))
    foreach(x -> _add_constraint!(m, imodel, x),MOI.get(imodel, MOI.ListOfConstraintIndices{SQF{T}, ET{T}}()))
    return nothing
end

function _add_conic_constraints!(opt::T, ip) where T
    foreach(fs -> MOI.add_constraint(opt, fs), values(_conic_socp(ip)))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Creates an initial node with initial box constraints and adds it to the stack.
"""
function create_initial_node!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    push!(m._stack, NodeBB{N,T}(ntuple(i -> _lower_bound(BranchVar, m, i), N),
                                ntuple(i -> _upper_bound(BranchVar, m, i), N),
                                -Inf, Inf, 1, 1, 0, BD_NONE))
    m._node_count = 1
    m._maximum_node_id += 1

    return nothing
end
function presolve_global!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    wp = _working_problem(m)
    ip = _input_problem(m)

    _add_decision_variables!(wp, ip)
    _add_linear_constraints!(wp, ip)
    _add_quadratic_constraints!(wp, ip)
    #_add_conic_constraints!(wp, ip)

    # set objective function
    m._working_problem._objective_sv = ip._objective_sv

    # add nonlinear constraints
    # the nonlinear evaluator loads with populated subexpressions which are then used
    # to asssess the linearity of subexpressions
    add_nonlinear_evaluator!(m)
    add_nonlinear_functions!(m)
    add_subexpression_buffers!(m)

    # converts a maximum problem to a minimum problem (internally) if necessary
    # this is placed after adding nonlinear functions as this prior routine
    # copies the nlp_block from the input_problem to the working problem
    _max_to_min!(m)
    reform_epigraph!(m)

    # labels variables to branch on
    label_branch_variables!(m)

    load_relaxed_problem!(m)
    create_initial_node!(m)

    branch_variable_num = m._branch_variable_num

    m._current_xref             = fill(0.0, branch_variable_num)
    m._candidate_xref           = fill(0.0, branch_variable_num)
    m._current_objective_xref   = fill(0.0, branch_variable_num)
    m._prior_objective_xref     = fill(0.0, branch_variable_num)
    m._lower_lvd                = fill(0.0, branch_variable_num)
    m._lower_uvd                = fill(0.0, branch_variable_num)

    # populate in full space until local MOI nlp solves support constraint deletion
    # uses input model for local nlp solves... may adjust this if a convincing reason
    # to use a reformulated upper problem presents itself
    m._lower_result_count_max  = 100
    m._lower_result_count      = 0
    m._lower_solution          = Vector{Float64}[]
    for j = 1:_lower_result_count_max
        push!(m._lower_solution, zeros(m._working_problem._variable_num))
    end

    m._cut_solution        = zeros(Float64, m._working_problem._variable_num)
    m._solution            = zeros(Float64, m._working_problem._variable_num)
    m._upper_solution      = zeros(Float64, m._working_problem._variable_num)
    m._upper_variables     = fill(VI(-1), m._working_problem._variable_num)

    # add storage for fbbt
    m._lower_fbbt_buffer   = zeros(Float64, m._working_problem._variable_num)
    m._upper_fbbt_buffer   = zeros(Float64, m._working_problem._variable_num)

    # add storage for obbt ( perform obbt on all relaxed variables, potentially)
    m._obbt_working_lower_index = fill(false, branch_variable_num)
    m._obbt_working_upper_index = fill(false, branch_variable_num)
    m._old_low_index            = fill(false, branch_variable_num)
    m._old_upp_index            = fill(false, branch_variable_num)
    m._new_low_index            = fill(false, branch_variable_num)
    m._new_upp_index            = fill(false, branch_variable_num)
    m._lower_indx_diff          = fill(false, branch_variable_num)
    m._upper_indx_diff          = fill(false, branch_variable_num)
    m._obbt_variable_num        = branch_variable_num

    # set subgradient refinement flag
    wp._relaxed_evaluator.is_post = m.subgrad_tighten
    wp._relaxed_evaluator.subgrad_tighten = m.subgrad_tighten
    wp._relaxed_evaluator.reverse_subgrad_tighten =  m.reverse_subgrad_tighten

    m._presolve_time = time() - m._parse_time

    return nothing
end
