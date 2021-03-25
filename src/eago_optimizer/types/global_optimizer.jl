struct FullVar end
struct BranchVar end

Base.@kwdef mutable struct GlobalOptimizer{N,T<:Real,S<:ExtensionType} <: MOI.AbstractOptimizer

        _ext_type::S = S()
        _stack::BinaryMinMaxHeap{NodeBB{N,T}} = BinaryMinMaxHeap{NodeBB{N,T}}()
        _current_node::NodeBB{N,T} = NodeBB{N,T}()

        _input_to_local_map::MOIU.IndexMap              = MOIU.IndexMap()
        _constraint_primal::Dict{CI,Union{T,Vector{T}}} = Dict{CI,T}()
        _constraint_offset::Vector{Int}                 = Int[]
        _constraint_index_num::Int = 0
        _constraint_row_num::Int = 0
        # loaded from _input_problem by TODO
        #_working_problem::ParsedProblem = ParsedProblem()

        _first_relax_point_set::Bool = false
        _current_xref::Vector{T} = T[]
        _candidate_xref::Vector{T} = T[]

        _use_prior_objective_xref::Bool = false
        _current_objective_xref::Vector{T} = T[]
        _prior_objective_xref::Vector{T} = T[]

        _user_branch_variables::Bool = false
        _fixed_variable::Vector{Bool} = Bool[]
        _branch_variable_num::Int = 0
        _branch_to_variable_map::Vector{Int} = Int[]
        _variable_to_branch_map::Vector{Int} = Int[]

        _lower_variable_bound::Vector{T} = T[]
        _upper_variable_bound::Vector{T} = T[]

        # all subproblem immutable subproblem status are set in global_solve in corresponding routines
        # in optimize_nonconvex.jl
        _preprocess_feasibility::Bool = true
        _preprocess_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
        _preprocess_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED

        _lower_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
        _lower_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
        _lower_feasibility::Bool = true
        _lower_objective_value::T = -Inf

        # set in TODO
        _lower_solution::Vector{T} = T[]
        _lower_lvd::Vector{T}      = T[]
        _lower_uvd::Vector{T}      = T[]

        _cut_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
        _cut_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
        _cut_objective_value::T = -Inf
        _cut_feasibility::Bool = true

        # set in TODO
        _cut_solution::Vector{T} = T[]

        _upper_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
        _upper_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
        _upper_feasibility::Bool = true
        _upper_objective_value::T = Inf

        # array is initialized to correct size in TODO, reset in single_nlp_solve! in optimize_convex.jl
        _upper_variables::Vector{VI} =  VI[]

        # set in TODO
        _upper_solution::Vector{T} = T[]

        _postprocess_feasibility::Bool = true

        _time_left::T = 1000.0
        _start_time::T = 0.0
        _run_time::T = 0.0
        _parse_time::T = 0.0
        _presolve_time::T = 0.0
        _last_preprocess_time::T = 0.0
        _last_lower_problem_time::T = 0.0
        _last_upper_problem_time::T = 0.0
        _last_postprocessing_time::T = 0.0

        # reset in initial_parse! in parse.jl
        _min_converged_value::T = Inf
        _global_lower_bound::T = -Inf
        _global_upper_bound::T = Inf
        _maximum_node_id::Int = 0
        _iteration_count::Int = 0
        _node_count::Int = 0

        # Storage for output, reset in initial_parse! in parse.jl
        _solution_value::T = 0.0
        _feasible_solution_found::Bool = false
        _first_solution_node::Int = -1
        _best_upper_value::T = Inf

        # Optimality-Based Bound Tightening (OBBT) Options
        # set in TODO
        _obbt_working_lower_index::Vector{Bool} = Bool[]
        _obbt_working_upper_index::Vector{Bool} = Bool[]
        _lower_indx_diff::Vector{Bool} = Bool[]
        _upper_indx_diff::Vector{Bool} = Bool[]
        _old_low_index::Vector{Bool} = Bool[]
        _old_upp_index::Vector{Bool} = Bool[]
        _new_low_index::Vector{Bool} = Bool[]
        _new_upp_index::Vector{Bool} = Bool[]
        _obbt_variables::Vector{VI} = VI[]
        _obbt_variable_num::Int = 0
        _obbt_performed_flag::Bool = false

        # Feasibility-Based Bound Tightening Options
        # set in set_constraint_propagation_fbbt in domain_reduction.jl
        _cp_improvement::T = 0.0
        _cp_evaluation_reverse::Bool = false

        _cut_iterations::Int = 0
        _cut_add_flag::Bool = false

        # Options for Repetition (If DBBT Performed Well)
        # set in within preprocess in optimize_nonconvex.jl
        _node_repetitions::Int = 0
        _initial_volume::T = 0.0
        _final_volume::T = 0.0

        # Log
        _log::Log = Log()

        # set in TODO
        _buffered_quadratic_ineq_ci::Vector{CI{SAF,LT}} = CI{SAF,LT}[]
        _buffered_quadratic_eq_ci::Vector{CI{SAF,LT}} = CI{SAF,LT}[]
        _buffered_nonlinear_ci::Vector{CI{SAF,LT}} = CI{SAF,LT}[]

        # initialized to empty in constructor (or via MOI.empty), filled in objective_cut in relax.jl
        # called by obbt in domain_reduction.jl, lower_problem, and add_cut in optimize_nonconvex.jl,
        # emptied in delete_objective_cuts! in relax.jl
        _objective_cut_ci_saf::Vector{CI{SAF,LT}} = CI{SAF,LT}[]

        # need to retreive primal _relaxed_variable_index
        # set in TODO
        #"Number of variables actively branched on in B&B routine (excludes linear and fixed)"
        _relaxed_variable_number::Int = 0
        _relaxed_variable_index::Vector{VI} = VI[]
        _relaxed_variable_eq::Vector{Tuple{CI{SV, ET}, Int}} = Tuple{CI{SV, ET}, Int}[]
        _relaxed_variable_lt::Vector{Tuple{CI{SV, LT}, Int}} = Tuple{CI{SV, LT}, Int}[]
        _relaxed_variable_gt::Vector{Tuple{CI{SV, GT}, Int}} = Tuple{CI{SV, GT}, Int}[]

        # set as user-input
        _branch_variables::Vector{Bool} = Bool[]

        _new_eval_constraint::Bool = false
        _new_eval_objective::Bool = false

        _node_to_sv_leq_ci::Vector{CI{SV,LT}} = CI{SV,LT}[]
        _node_to_sv_geq_ci::Vector{CI{SV,GT}} = CI{SV,GT}[]

        #"Set to true if a nonlinear evaluator was created (NLconstraint or NLobjective specified)"
        _nonlinear_evaluator_created::Bool = false

        #_relaxed_evaluator::Evaluator = Evaluator{1,NS}()
        #_relaxed_constraint_bounds::Vector{MOI.NLPBoundsPair} = Vector{MOI.NLPBoundsPair}[]
end
GlobalOptimizer(::Type{T}) where T<:AbstractFloat = GlobalOptimizer{1,T,DefaultExt}()

function MOI.is_empty(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    empty_opt = GlobalOptimizer{N,T,S}()
    is_empty_flag = true
    for f in fieldnames(GlobalOptimizer)
        if f == :_stack
            is_empty_flag &= isempty(getfield(m, f))
        elseif f == :_current_node
            is_empty_flag &= getfield(m, f) == NodeBB{N,T}()
        elseif f == :_log
        else
            if getfield(empty_opt, f) != getfield(m, f)
                is_empty_flag = false
                break
            end
        end
    end
    return is_empty_flag
end

#=
TODO: StructArray access is faster here for most...
=#

Base.@propagate_inbounds function _is_unfixed_integer(::FullVar, m::GlobalOptimizer, i)
end
Base.@propagate_inbounds function _is_unfixed_integer(::BranchVar, m::GlobalOptimizer, i)
end

for (f,arr) in ((:_lower_bound, :_lower_variable_bound),
                (:_upper_bound, :_upper_variable_bound),
                (:_lower_solution, :_lower_solution),)
    @eval function ($f)(::FullVar, m::GlobalOptimizer{N,T,S}) where {N,T<:Real,S<:ExtensionType}
        m.$arr
    end
    @eval Base.@propagate_inbounds function ($f)(::FullVar, m::GlobalOptimizer{N,T,S}, i) where {N,T<:Real,S<:ExtensionType}
        m.$arr[i]
    end
    @eval Base.@propagate_inbounds function ($f)(::BranchVar, m::GlobalOptimizer{N,T,S}) where {N,T<:Real,S<:ExtensionType}
        m.$arr[m._branch_variable]
    end
    @eval Base.@propagate_inbounds function ($f)(::BranchVar, m::GlobalOptimizer{N,T,S}, i) where {N,T<:Real,S<:ExtensionType}
        m.$arr[m._branch_variable[i]]
    end
end

Base.@propagate_inbounds function _set_lower_bound!(::BranchVar, m::GlobalOptimizer{N,T,S}, i, v) where {N,T<:Real,S<:ExtensionType}
    m._lower_variable_bound[m._branch_variable[i]] = v
end
Base.@propagate_inbounds function _set_lower_bound!(::FullVar, m::GlobalOptimizer{N,T,S}, i, v) where {N,T<:Real,S<:ExtensionType}
    m._lower_variable_bound[i] = v
end
Base.@propagate_inbounds function _set_upper_bound!(::BranchVar, m::GlobalOptimizer{N,T,S}, i, v) where {N,T<:Real,S<:ExtensionType}
    m._upper_variable_bound[m._branch_variable[i]] = v
end
Base.@propagate_inbounds function _set_upper_bound!(::FullVar, m::GlobalOptimizer{N,T,S}, i, v) where {N,T<:Real,S<:ExtensionType}
    m._upper_variable_bound[i] = v
end


function _mid(::FullVar, m::GlobalOptimizer{N,T,S}) where {N,T<:Real,S<:ExtensionType}
    0.5*(_upper_bound(FullVar, m) - _lower_bound(FullVar, m))
end
Base.@propagate_inbounds function _mid(::FullVar, m::GlobalOptimizer{N,T,S}, i) where {N,T<:Real,S<:ExtensionType}
    0.5*(_upper_bound(FullVar, m, i) - _lower_bound(FullVar, m, i))
end
function _mid(::BranchVar, m::GlobalOptimizer{N,T,S}) where {N,T<:Real,S<:ExtensionType}
end
Base.@propagate_inbounds function _mid(::BranchVar, m::GlobalOptimizer{N,T,S}, i) where {N,T<:Real,S<:ExtensionType}
end
