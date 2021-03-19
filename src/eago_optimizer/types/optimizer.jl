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
Base.@kwdef mutable struct InputModel
    _input_model::MOIU.Model{Float64} = MOIU.Model{Float64}()
    _nlp_data::Union{MOI.NLPBlockData, Nothing} = nothing
end

@enum(ProblemType, UNCLASSIFIED, LP, MILP, SOCP, MISOCP, SDP, DIFF_CVX, MICVX, MINCVX)

#=
Holds outputs in the Optimizer object, passing parameters EAGO specific parameters
to _options, the GlobalOptimizer is remade to ensure type stability of subroutines
prior to initializing a global solve.
=#
export Optimizer
Base.@kwdef mutable struct Optimizer <: MOI.AbstractOptimizer

    _model::InputProblem             = InputModel()
    _options::GlobalOptimizerOptions = GlobalOptimizerOptions()
    _solver::GlobalOptimizer         = GlobalOptimizer()

    _problem_type::ProblemType = UNCLASSIFIED

    _termination_status_code::MOI.TerminationStatusCode = MOI.OPTIMIZE_NOT_CALLED
    _primal_status_code::MOI.ResultStatusCode           = MOI.OTHER_RESULT_STATUS
    _dual_status_code::MOI.ResultStatusCode             = MOI.OTHER_RESULT_STATUS

    _objective_value::Float64        = Inf
    _objective_bound::Float64        = -Inf
    _iteration_count::Int64          = 0
    _node_count::Int64               = 0

    _is_silent::Bool     = false

    _start_time::Float64   = 0.0
    _time_limit::Float64 = Inf
    _time_left::Float64  = 1000.0
    _run_time::Float64   = 0.0
    _parse_time::Float64   = 0.0

    _solution::Vector{Float64}                = Float64[]
    _constraint_primal::Vector{Float64}       = Float64[]
end


for f in fieldnames(GlobalOptimizerOptions)
    @eval $(Symbol("_"*String(F)))(x::Optimizer) = x._options.$F
end

function _set_cons_primal!(m::Optimizer, ci::CI{F,S}, v::Float64) where{F <: MOI.AbstractFunction,
                                                                        S <: MOI.AbstractScalarSet}
    m._constraint_primal[ci] = v
    return
end

function _set_cons_primal!(m::Optimizer, ci::CI{F,S}, v::Vector{Float64}) where{F <: MOI.AbstractFunction,
                                                                                S <: MOI.AbstractVectorSet}
    m._constraint_primal[ci] = v
    return
end

@inline function _get_variable_lo(::typeof(:branch), m::Optimizer, i)
    return _lower_variable_bound(m._current_node, m._sol_to_branch_map[i])
end

@inline function _get_variable_hi(::typeof(:branch), m::Optimizer, i)
    return _upper_variable_bound(m._current_node, m._sol_to_branch_map[i])
end
