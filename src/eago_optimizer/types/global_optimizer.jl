Base.@kwdef mutable struct GlobalOptimizer <: MOI.AbstractOptimizer

        _input_to_local_map::MOIU.IndexMap                            = MOIU.IndexMap()
        _constraint_primal::Dict{CI,Union{Float64,Vector{Float64}}}   = Dict{CI,Float64}()
        _constraint_offset::Vector{Int}             = Int[]
        _constraint_index_num::Int = 0
        _constraint_row_num::Int = 0
        # loaded from _input_problem by TODO
        #_working_problem::ParsedProblem = ParsedProblem()

        _stack::BinaryMinMaxHeap{NodeBB} = BinaryMinMaxHeap{NodeBB}()

        # set in node_selection!
        _current_node::NodeBB = NodeBB()

        _first_relax_point_set::Bool = false
        _current_xref::Vector{Float64} = Float64[]
        _candidate_xref::Vector{Float64} = Float64[]

        _use_prior_objective_xref::Bool = false
        _current_objective_xref::Vector{Float64} = Float64[]
        _prior_objective_xref::Vector{Float64} = Float64[]

        # set in label_branch_variables! and label_fixed_variables! respectively in parse.jl
        _user_branch_variables::Bool = false
        _fixed_variable::Vector{Bool} = Bool[]
        _branch_variable_num::Int = 0
        _branch_to_sol_map::Vector{Int} = Int[]
        _sol_to_branch_map::Vector{Int} = Int[]

        # all subproblem immutable subproblem status are set in global_solve in corresponding routines
        # in optimize_nonconvex.jl
        _preprocess_feasibility::Bool = true
        _preprocess_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
        _preprocess_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED

        _lower_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
        _lower_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
        _lower_feasibility::Bool = true
        _lower_objective_value::Float64 = -Inf

        # set in TODO
        _lower_solution::Vector{Float64} = Float64[]
        _lower_lvd::Vector{Float64} = Float64[]
        _lower_uvd::Vector{Float64} = Float64[]

        _cut_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
        _cut_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
        _cut_objective_value::Float64 = -Inf
        _cut_feasibility::Bool = true

        # set in TODO
        _cut_solution::Vector{Float64} = Float64[]

        _upper_result_status::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS
        _upper_termination_status::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
        _upper_feasibility::Bool = true
        _upper_objective_value::Float64 = Inf

        # array is initialized to correct size in TODO, reset in single_nlp_solve! in optimize_convex.jl
        _upper_variables::Vector{VI} =  VI[]

        # set in TODO
        _upper_solution::Vector{Float64} = Float64[]

        _postprocess_feasibility::Bool = true

        _time_left::Float64 = 1000.0
        _start_time::Float64 = 0.0
        _run_time::Float64 = 0.0
        _parse_time::Float64 = 0.0
        _presolve_time::Float64 = 0.0
        _last_preprocess_time::Float64 = 0.0
        _last_lower_problem_time::Float64 = 0.0
        _last_upper_problem_time::Float64 = 0.0
        _last_postprocessing_time::Float64 = 0.0

        # reset in initial_parse! in parse.jl
        _min_converged_value::Float64 = Inf
        _global_lower_bound::Float64 = -Inf
        _global_upper_bound::Float64 = Inf
        _maximum_node_id::Int64 = 0
        _iteration_count::Int64 = 0
        _node_count::Int64 = 0

        # Storage for output, reset in initial_parse! in parse.jl
        _solution_value::Float64 = 0.0
        _feasible_solution_found::Bool = false
        _first_solution_node::Int64 = -1
        _best_upper_value::Float64 = Inf

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

        # Buffers for fbbt, set in presolve, used in preprocess
        _lower_fbbt_buffer::Vector{Float64} = Float64[]
        _upper_fbbt_buffer::Vector{Float64} = Float64[]

        # Feasibility-Based Bound Tightening Options
        # set in set_constraint_propagation_fbbt in domain_reduction.jl
        _cp_improvement::Float64 = 0.0
        _cp_evaluation_reverse::Bool = false

        _cut_iterations::Int64 = 0
        _cut_add_flag::Bool = false

        # Options for Repetition (If DBBT Performed Well)
        # set in within preprocess in optimize_nonconvex.jl
        _node_repetitions::Int64 = 0
        _initial_volume::Float64 = 0.0
        _final_volume::Float64 = 0.0

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

function MOI.is_empty(m::GlobalOptimizer)
    empty_opt = GlobalOptimizer()
    is_empty_flag = true
    for f in fieldnames(GlobalOptimizer)
        if f == :_stack
            is_empty_flag &= isempty(getfield(m, f))
        elseif f == :_current_node
            is_empty_flag &= getfield(m, f) != NodeBB()
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
