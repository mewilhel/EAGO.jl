export EAGOParameters
"""
$(TYPEDEF)

Storage for parameters that do not change during a global solve.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct EAGOParameters

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
end