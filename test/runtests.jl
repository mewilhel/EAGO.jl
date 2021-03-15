#!/usr/bin/env julia

using Test
using EAGO, JuMP, MathOptInterface, Ipopt, ForwardDiff
using IntervalArithmetic, SpecialFunctions
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
using ForwardDiff: Dual, Partials

unit_excludes = String[
    "number_threads",                # EAGO won't support number of threads in near future
    "solve_qcp_edge_cases",
    "raw_status_string",             # TODO: ADD internal status states to EAGO
    "solve_qp_zero_offdiag",
    #"variablenames",                # SEEMS LIKE A FALSE PASS
    "solve_integer_edge_cases",
    "solve_zero_one_with_bounds_3",
    "solve_zero_one_with_bounds_2",
    "solve_qp_edge_cases",
    "solve_result_index",
    "solve_objbound_edge_cases",
]

contlinear_excludes = String[
    "linear14",
    "linear4",
    "linear1",
    "partial_start",
    "linear13"
]

intlinear_excludes = String[
    #=
    "int1",
    "int2",
    "int3",
    "indicator1",
    "indicator2",
    "indicator3",
    "indicator4",
    "semiinttest",
    "knapsack",
    "semiconttest"
    =#
]

contconic_excludes = String[
    "normone",
    "norminf",
    "sdp",
    "normnuc",
    "dualexp",
    "dualpow",
    "exp",
    "lin",
    "logdet",
    "soc",
    "normspec",
    "rootdet",
    "relentr",
    "rsoc",
    "pow",
    "geomean"
]

contquadratic_excludes = String[
    "ncqcp",
    "qp",
    "socp",
    "qcp",
]

nlp_excludes = String[
    "feasibility_sense_with_objective_and_hessian",
    "feasibility_sense_with_no_objective_and_no_hessian",
    "feasibility_sense_with_objective_and_no_hessian",
    "hs071_no_hessian",
    "feasibility_sense_with_no_objective_and_with_hessian",
    "nlp_objective_and_moi_objective",
    "hs071"
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

#include("branch_bound.jl")
#include("domain_reduction.jl")
#include("optimizer.jl")
#include("script_optimizer.jl")
#include("semiinfinite.jl")

default_options = (:verbosity => 0,)
test_time = @elapsed test_moi(Float64; default_options...)
#@printf("%8.2e seconds\n", test_time)
