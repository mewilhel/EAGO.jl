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
    socp_optimizer::MOI.AbstractOptimizer               = COSMO.Optimizer()
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

    # Default optimizer parameters
    # Presolving options
    "Should EAGO attempt to remove type assert issues for user-defined functions (default = false)"
    presolve_scrubber_flag::Bool = false
    "Create and use DAG representations of user-defined function (default = false)."
    presolve_to_JuMP_flag::Bool = false
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Apply the epigraph reformulation
    to the problem (default = false)."
    presolve_epigraph_flag::Bool = false
    "Rerranges the DAG using registered transformations (default = false)"
    presolve_flatten_flag::Bool = false

    # Conic reformulations
    "Attempt to bridge convex constraint to second order cone"
    conic_convert_quadratic::Bool = false

    # Iteration logging options
    "Turns logging on records global bounds, node count and run time. Additional
     options are available for recording information specific to subproblems (default = false)."
    log_on::Bool = false
    "Turns on logging of times and feasibility of subproblems (default = false)"
    log_subproblem_info::Bool = false
    "Log data every `log_interval` iterations (default = 1)."
    log_interval::Int64 = 1

    # Optimizer display options
    "The amount of information that should be printed to console while solving
    values range from 0 - 4: 0 is silent, 1 shows iteration summary statistics
    only, 2-4 show varying degrees of details about calculations within each
    iteration (default = 1)."
    verbosity::Int64 = 1
    "Display summary of iteration to console every `output_iterations` (default = 10)"
    output_iterations::Int64 = 1000
    "Display header for summary to console every `output_iterations` (default = 100)"
    header_iterations::Int64 = 10000

    # Node branching options
    "Convex coefficient used to select branch point. Branch point is given by
    `branch_cvx_factor*xmid + (1-branch_cvx_factor)*xsol` (default = 0.25)"
    branch_cvx_factor::Float64 = 0.25
    "Minimum distance from bound to have branch point normalized by width of
    dimension to branch on (default = 0.15)"
    branch_offset::Float64 = 0.15
    "Variables to branch on (default is all nonlinear)."
    branch_variable::Vector{Bool} = Bool[]
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Number of times repeat node
    processing priorto branching (default = 4)."
    branch_max_repetitions::Int64 = 4
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Volume ratio tolerance required
    to repeat processing the current node (default = 0.9)"
    branch_repetition_tol::Float64 = 0.9

    # Termination limits
    "Maximum number of nodes (default = 1E-7)"
    node_limit::Int64 = 1*10^7
    "Maximum CPU time in seconds (default = 1000)"
    time_limit::Float64 = 1000.0
    "Maximum number of iterations (default 3E6)"
    iteration_limit::Int64 = 3*10^6
    "Absolute tolerance for termination (default = 1E-3)"
    absolute_tolerance::Float64 = 1E-3
    "Relative tolerance for termination (default = 1E-3)"
    relative_tolerance::Float64 = 1E-3
    "Absolute constraint feasibility tolerance"
    absolute_constraint_feas_tolerance::Float64 = 1E-6
    "Perform only a local solve of the problem (default = false)."
    local_solve_only::Bool = false
    "[TO BE REMOVED] Flag stops B&B loop if feasible point found."
    feasible_local_continue::Bool = false

    # Options for constraint propagation
    "Depth in B&B tree above which constraint propagation should be disabled (default = 1000)"
    cp_depth::Int64 = 20
    "Number of times to repeat forward-reverse pass routine (default = 3)"
    cp_repetitions::Int64 = 3
    "Disable constraint propagation if the ratio of new node volume to beginning node volume exceeds
    this number (default = 0.99)"
    cp_tolerance::Float64 = 0.99
    "Use only valid interval bounds during constraint propagation (default = false)"
    cp_interval_only::Bool = false

    # obbt options
    "Depth in B&B tree above which OBBT should be disabled (default = 6)"
    obbt_depth::Int64 = 4
    "Number of repetitions of OBBT to perform in preprocessing (default = 3)"
    obbt_repetitions::Int64 = 20
    "Turn aggresive OBBT on (default = false)"
    obbt_aggressive_on::Bool = true
    "Maximum iteration to perform aggresive OBBT (default = 2)"
    obbt_aggressive_max_iteration::Int64 = 2
    "Minimum dimension to perform aggresive OBBT (default = 2)"
    obbt_aggressive_min_dimension::Int64 = 2
    "Tolerance to consider bounds equal (default = 1E-9)"
    obbt_tolerance::Float64 = 1E-9

    # Options for linear bound tightening
    "Depth in B&B tree above which linear FBBT should be disabled (default = 1000)"
    fbbt_lp_depth::Int64  = 1000
    "Number of repetitions of linear FBBT to perform in preprocessing (default = 3)"
    fbbt_lp_repetitions::Int64  = 3

    # Options for quadratic bound tightening
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Depth in B&B tree above which univariate quadratic FBBT should be disabled (default = -1)"
    quad_uni_depth::Int64 = -1
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Number of repetitions of univariate quadratic FBBT to perform in preprocessing (default = 2)"
    quad_uni_repetitions::Int64 = 2
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Depth in B&B tree above which bivariate
    quadratic FBBT should be disabled (default = -1)"
    quad_bi_depth::Int64 = -1
    "[FUTURE FEATURE, NOT CURRENTLY IMPLEMENTED] Number of repetitions of bivariate quadratic FBBT to perform in preprocessing (default = 2)."
    quad_bi_repetitions::Int64 = 2

    # Duality-based bound tightening (DBBT) options
    "Depth in B&B tree above which duality-based bound tightening should be disabled (default = 1E10)"
    dbbt_depth::Int64 = 10^10
    "New bound is considered equal to the prior bound if within dbbt_tolerance (default = 1E-9)."
    dbbt_tolerance::Float64 = 1E-8

    # Subgradient tightening flag
    "Relax Tag used to specify type of McCormick operator"
    relax_tag::RelaxTag = NS()
    "Perform tightening of interval bounds using subgradients at each factor in
    each nonlinear tape during a forward pass (default = true)."
    subgrad_tighten::Bool = true
    "Perform tightening of interval bounds using subgradients at each factor in
    each nonlinear tape during a reverse pass (default = false)."
    reverse_subgrad_tighten::Bool = false
    "Outer round computed subgradient bounds by this amount"
    subgrad_tol::Float64 = 1E-10

    # Tolerance to add cuts and max number of cuts
    "Minimum number of cuts at each node to attempt (unsafe cuts not necessarily added)"
    cut_min_iterations::Int64 = 1
    "Maximum number of cuts at each node to attempt"
    cut_max_iterations::Int64 = 3
    "Convex coefficient used to select point for new added cuts. Branch point is
    given by `(1-cut_cvx)*xmid + cut_cvx*xsol` (default = 0.9)."
    cut_cvx::Float64 = 0.9
    "Add cut if the L1 distance from the prior cutting point to the new cutting
    point normalized by the box volume is greater than the tolerance (default = 0.05)."
    cut_tolerance::Float64 = 0.05
    "Adds an objective cut to the relaxed problem (default = true)."
    objective_cut_on::Bool = true

    "Use tolerances to determine safe cuts in a Khajavirad 2018 manner"
    cut_safe_on::Bool = true
    "Lower tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_l::Float64 = 1E-8
    "Upper tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_u::Float64 = 1E8
    "Constant tolerance for safe-lp cut, Khajavirad 2018"
    cut_safe_b::Float64 = 1E9

    "Solve upper problem for every node with depth less than `upper_bounding_depth`
    and with a probabilityof (1/2)^(depth-upper_bounding_depth) otherwise (default = 6)"
    upper_bounding_depth::Int64 = 6

    # handling for domain violations
    "Amount about a domain violation to ignore when propagating bounds."
    domain_violation_guard_on::Bool = false
    "Amount about a domain violation to ignore when propagating bounds."
    domain_violation_ϵ::Float64 = 1E-9

    # set by MOI manipulations
    _input_problem::MOIU.Model{Float64} = MOIU.Model{Float64}()
    _nlp_data::Union{MOI.NLPBlockData, Nothing} = nothing
    _constraint_offset::Vector{Int}      = Int[]
    _constraint_index_num::Int = 0
    _constraint_row_num::Int = 0

    _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _result_status_code::MOI.ResultStatusCode = MOI.OTHER_RESULT_STATUS

    _solution::Vector{Float64}                = Float64[]
    _primal_constraint_value::Vector{Float64} = Float64[]

    # loaded from _input_problem by TODO
    _working_problem::ParsedProblem = ParsedProblem()

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
    _objective_value::Float64 = Inf
    _objective_bound::Float64 = -Inf
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

_input_problem(m::Optimizer) = m._input_problem
_working_problem(m::Optimizer) = m._working_problem

function _set_cons_primal!(m::Optimizer, ci::CI{F,S}, v::Float64) where{F <: MOI.AbstractFunction,
                                                                        S <: MOI.AbstractScalarSet}
    i = ci.value
    os = m._constraint_offset
    m._primal_constraint_value[os[i]] = v
    return
end

function _set_cons_primal!(m::Optimizer, ci::CI{F,S}, v::Vector{Float64}) where{F <: MOI.AbstractFunction,
                                                                                S <: MOI.AbstractVectorSet}
    i = ci.value
    os = m._constraint_offset
    m._primal_constraint_value[(os[i] + 1):os[i + 1]] .= v
    return
end
