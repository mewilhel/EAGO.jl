#!/usr/bin/env julia

using Test
using EAGO, JuMP, MathOptInterface, Ipopt, ForwardDiff
using IntervalArithmetic, SpecialFunctions
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
using ForwardDiff: Dual, Partials


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

    @testset "linear tests" begin
        MOIT.contlineartest(optimizer, config)
    end
end

#include("branch_bound.jl")
#include("domain_reduction.jl")
#include("optimizer.jl")
#include("script_optimizer.jl")
#include("semiinfinite.jl")

default_options = (:verbosity => 0,)
test_time = test_moi(Float64; default_options...)
@printf("%8.2e seconds\n", test_time)
