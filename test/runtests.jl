#!/usr/bin/env julia

using Test, Printf, EAGO, MathOptInterface #, MINLPTests

const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

include(joinpath(@__DIR__, "moit_tests.jl"))
#include(joinpath(@__DIR__, "minlp_tests.jl"))

#include(joinpath(@__DIR__, "branch_bound.jl"))
#include(joinpath(@__DIR__, "domain_reduction.jl"))
#include(joinpath(@__DIR__, "semiinfinite_tests.jl"))

#, JuMP, Ipopt, ForwardDiff
#using IntervalArithmetic, SpecialFunctions
#using ForwardDiff: Dual, Partials
#include("script_optimizer.jl")
