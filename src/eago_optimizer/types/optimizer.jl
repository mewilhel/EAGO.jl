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

# Holds a copy of the problem passed to EAGO
# Needs separate storage for the MOI compilant model and the nlp_data
Base.@kwdef mutable struct InputModel{T<:Real}
    _input_model::MOIU.Model{T} = MOIU.Model{T}()
    _nlp_data::Union{MOI.NLPBlockData, Nothing} = nothing
end

MOI.is_empty(m::InputModel) = (m._nlp_data == nothing) && MOI.is_empty(m._input_model)

@enum(ProblemType, UNCLASSIFIED, LP, MILP, SOCP, MISOCP, SDP, DIFF_CVX, MICVX, MINCVX)

#=
Holds outputs in the Optimizer object, passing parameters EAGO specific parameters
to _options, the GlobalOptimizer is remade to ensure type stability of subroutines
prior to initializing a global solve.
=#
export Optimizer
mutable struct Optimizer{T<:AbstractFloat} <: MOI.AbstractOptimizer
    _model::InputModel{T}
    _options::GlobalOptimizerOptions{T}
    _solver::GlobalOptimizer
    _problem_type::ProblemType

    _termination_status_code::MOI.TerminationStatusCode
    _primal_status_code::MOI.ResultStatusCode

    _objective_value::T
    _objective_bound::T
    _relative_gap::T
    _iteration_count::Int
    _node_count::Int

    _is_silent::Bool

    _start_time::T
    _time_limit::T
    _time_left::T
    _run_time::T
    _parse_time::T

    _solution::Vector{T}
    _constraint_primal::Dict{CI,Union{Vector{T},T}}
    _input_to_solution_map
end
function Optimizer{T}(; kwargs...) where T<:Real
    model   = InputModel{T}()
    options = GlobalOptimizerOptions{T}(;kwargs...)
    solver  = GlobalOptimizer(T)
    problem_type = UNCLASSIFIED

    termination_status_code = MOI.OPTIMIZE_NOT_CALLED
    primal_status_code      = MOI.OTHER_RESULT_STATUS

    objective_value = Inf
    objective_bound = -Inf
    relative_gap = Inf
    iteration_count = 0
    node_count      = 0

    is_silent       = false

    start_time   = 0.0
    time_limit   = Inf
    time_left    = 1000.0
    run_time     = 0.0
    parse_time   = 0.0

    solution          = T[]
    constraint_primal = Dict{CI,T}()
    input_to_solution_map = nothing

    return Optimizer{T}(model,
                        options,
                        solver,
                        problem_type,
                        termination_status_code,
                        primal_status_code,
                        objective_value,
                        objective_bound,
                        relative_gap,
                        iteration_count,
                        node_count,
                        is_silent,
                        start_time,
                        time_limit,
                        time_left,
                        run_time,
                        parse_time,
                        solution,
                        constraint_primal,
                        input_to_solution_map)
end
Optimizer(; kwargs...) = Optimizer{Float64}(; kwargs...)

for f in fieldnames(GlobalOptimizerOptions)
    @eval @inline $(Symbol("_"*String(f)))(x::Optimizer) = x._options.$f
end

_input_model(m::Optimizer{T}) where T = m._model._input_model
_input_nlp_data(m::Optimizer) = m._model._nlp_data

function _set_cons_primal!(m::Optimizer{T}, ci::CI{F,S}, v::T) where {F <: MOI.AbstractFunction,
                                                                     S <: MOI.AbstractScalarSet,
                                                                     T <: AbstractFloat}
    m._constraint_primal[ci] = v
    return
end

function _set_cons_primal!(m::Optimizer{T}, ci::CI{F,S}, v::Vector{T}) where {F <: MOI.AbstractFunction,
                                                                              S <: MOI.AbstractVectorSet,
                                                                              T <: AbstractFloat}
    m._constraint_primal[ci] = v
    return
end
