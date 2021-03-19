
const unit_excludes = String[
    "number_threads",                 # EAGO won't support number of threads in near future
    "raw_status_string",              # TODO: ADD internal status states to EAGO
    "variablenames",                  # SEEMS LIKE A FALSE PASS

    "solve_qcp_edge_cases",           # Need to update nonconvex...
    "solve_qp_zero_offdiag",
     "solve_qp_edge_cases",

    "solve_integer_edge_cases",       # Integer valued variable not currently supported.
    "solve_objbound_edge_cases",

    "solve_zero_one_with_bounds_3",   # GLPK has a non-standard return code
    "solve_result_index",             # TODO: Should throw error when querying for multiple results... (expected behavior?)

    #"getvariable",
    #"solve_zero_one_with_bounds_2",
    #"getconstraint",
    #"solve_with_upperbound",
    #"solve_single_variable_dual_min",
    #"solve_single_variable_dual_max",
    #"solve_with_lowerbound",
    #"solve_zero_one_with_bounds_1"
]

const contlinear_excludes = String[
    "partial_start",                 # EAGO doesn't support VariablePrimalStart
    #"linear14",
    #"linear1"
]

const intlinear_excludes = String[
    "indicator1",    # doesn't currently support indicator sets
    "indicator2",
    "indicator3",
    "indicator4",

    "int2",
    "int3",

    "semiinttest",   # doesn't currently support SV in SemiInteger
    #"semiconttest"
]

const contconic_excludes = String[
    "sdp",
    "normnuc",
    "dualexp",
    "dualpow",
    "exp",
    "logdet",
    #"soc",
    "normspec",
    "rootdet",
    "relentr",
    "rsoc",
    "pow",
    "lin",
    "geomean"
]

const contquadratic_excludes = String[
    "ncqcp",
    "qp",
    "socp",
    "qcp",
]

const nlp_excludes = String[
    "feasibility_sense_with_objective_and_hessian",
    "feasibility_sense_with_no_objective_and_no_hessian",
    "feasibility_sense_with_objective_and_no_hessian",
    "hs071_no_hessian",
    "feasibility_sense_with_no_objective_and_with_hessian",
    "nlp_objective_and_moi_objective",
    "hs071",
    "hs071_hessian_vector_product_test"
]

function test_moi(T::Type{<:Real}; solver_options...)

    optimizer = MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{T}()), EAGO.Optimizer(; solver_options...))

    tol = 2sqrt(sqrt(eps(T)))
    config = MOIT.TestConfig{T}(
    atol = tol,
    rtol = tol,
    solve = true,
    query = true,
    modify_lhs = false,
    duals = false,
    infeas_certificates = false,
    )

    @testset "unit tests" begin
        MOIT.unittest(MOIB.full_bridge_optimizer(optimizer, T), config, unit_excludes)
    end

    @testset "continuous linear tests" begin
        MOIT.contlineartest(MOIB.full_bridge_optimizer(optimizer, T), config, contlinear_excludes)
    end

    @testset "mixed-integer linear tests" begin
        MOIT.intlineartest(MOIB.full_bridge_optimizer(optimizer, T), config, intlinear_excludes)
    end

    @testset "continuous conic tests" begin
        MOIT.contconictest(MOIB.full_bridge_optimizer(optimizer, T), config, contconic_excludes)
    end

    @testset "continuous quadratic tests" begin
        MOIT.contquadratictest(MOIB.full_bridge_optimizer(optimizer, T), config, contquadratic_excludes)
    end

    @testset "nlp" begin
        MOIT.nlptest(MOIB.full_bridge_optimizer(optimizer, T), config, nlp_excludes)
    end
end

default_options = (:verbosity => 0,)
test_moi(Float64; default_options...)
