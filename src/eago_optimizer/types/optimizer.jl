# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimizer.jl
# Defines optimizer structure used by EAGO. Namely, ObjectiveType, ProblemType
# EAGOParameters, InputProblem, ParsedProblem, and Optimizer.
#############################################################################

function default_nlp_solver()

    upper_optimizer = Ipopt.Optimizer()

    MOI.set(upper_optimizer, MOI.RawParameter("max_iter"),3000)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_tol"), 1E30)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_iter"), 300)
    MOI.set(upper_optimizer, MOI.RawParameter("constr_viol_tol"), 0.000001)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_compl_inf_tol"), 0.000001)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_dual_inf_tol"), 1.0)
    MOI.set(upper_optimizer, MOI.RawParameter("acceptable_constr_viol_tol"), 0.000001)
    MOI.set(upper_optimizer, MOI.RawParameter("print_level"), 0)

    return upper_optimizer
end

export Optimizer
"""
$(TYPEDEF)

The main optimizer object used by EAGO to solve problems during the optimization
routine. The following commonly used options are described below and can be set
via keyword arguments in the JuMP/MOI model. The raw parameter interface however
is likely preferable. The Optimizer is organized in the following manner. Parameters
which are expected to be constant over the entire solve are stored in
`_parameters::EAGOParameters` field. User-facing keywords not in EAGOParameters field:
- `relaxed_optimizer::MOI.AbstractOptimizer`: An instance of the optimizer used to solve the relaxed subproblems (default = GLPK.Optimizer())
- `obbt_variable_values::Vector{Bool}`: Variables to perform OBBT on (default: all variables in nonlinear expressions).
- `upper_optimizer::MOI.AbstractOptimizer`: Optimizer used to solve upper bounding problems. (default = Ipopt.Optimizer)
- `enable_optimize_hook::Bool`: Specifies that the optimize_hook! function should be called rather than throw the problem to the standard B&B routine (default = false).
- `ext::Dict{Symbol, Any}`: Holds additional storage needed for constructing extensions to EAGO (default = Dict{Symbol,Any}).
- `ext_type::ExtensionType`: Holds an instance of a subtype of `EAGO.ExtensionType` used to define new custom subroutines (default = DefaultExt()).
"""
Base.@kwdef mutable struct Optimizer <: MOI.AbstractOptimizer

    # Subsolver Options
    "Subsolver used to solve linear programs"
    lp_optimizer::MOI.AbstractOptimizer                 = GLPK.Optimizer()
    "Subsolver used to solve mixed-integer linear programs"
    mip_optimizer::MOI.AbstractOptimizer                = GLPK.Optimizer()
    "Subsolver used to solve second-order conic programs"
    second_order_cone_optimizer::MOI.AbstractOptimizer  = COSMO.Optimizer()
    "Semidefinite Programming Optimizer"
    semidefinite_optimizer::MOI.AbstractOptimizer       = COSMO.Optimizer()
    "Subsolver used to locally solve nonlinear programs"
    nlp_optimizer::MOI.AbstractOptimizer                = Ipopt.Optimizer()
    "Subsolver used to locally solve mixed-integer nonlinear programs "
    minlp_optimizer::MOI.AbstractOptimizer              = Pavito.Optimizer()

    # Options for optimality-based bound tightening
    # set as a user-specified option
    relaxed_optimizer::MOI.AbstractOptimizer = GLPK.Optimizer()

    # Upper bounding options (set as a user-specified option)
    upper_optimizer::MOI.AbstractOptimizer = default_nlp_solver()

    # set as a user-specified option (if empty set to all nonlinear by TODO in TODO)
    obbt_variable_values::Vector{Bool} = Bool[]

    # Extensions (set as user-specified option)
    enable_optimize_hook::Bool = false
    ext::Dict{Symbol, Any} = Dict{Symbol,Any}()
    ext_type::ExtensionType = DefaultExt()

    # set as user-specified option
    _parameters::EAGOParameters = EAGOParameters()

    # set by MOI manipulations (see Input problem structure)
    _input_problem::InputProblem = InputProblem()

    # loaded from _input_problem by TODO
    _working_problem::ParsedProblem = ParsedProblem()

    _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _result_status_code::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS

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
    _branch_variable_count::Int = 0
    _branch_to_sol_map::Vector{Int} = Int[]
    _sol_to_branch_map::Vector{Int} = Int[]

    _continuous_solution::Vector{Float64} = Float64[]

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

    # set to time limit in initial_parse! in parse.jl, decremented throughout global_solve in optimize_nonconvex.jl
    _time_left::Float64 = 1000.0

    # set constructor reset on empty! and  to zero in initial parse! in parse.jl
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
    _objective_value::Float64 = -Inf
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
    _obbt_variable_count::Int = 0
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

    # set initially in TODO, reset in objective_cut in relax.jl
    _objective_cut_ci_sv::CI{SV,LT} = CI{SV,LT}(-1)

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
