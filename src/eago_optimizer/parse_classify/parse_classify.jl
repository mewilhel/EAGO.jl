# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# Defines functions used to parse the input optimization problem into
# a solvable form including routines used to classify input problems as a
# LP, SOCP, MILP, MISOCP, and convex problem types.
#############################################################################

#include(joinpath(@__DIR__, "parse_bridge_nl", "parse_bridge_nl.jl"))
include(joinpath(@__DIR__, "label_branch.jl"))

const INTEGER_FUNC_SET{T <: Real} = Union{Tuple{MOI.SingleVariable, MOI.ZeroOne},
                                          Tuple{MOI.SingleVariable, MOI.Integer},
                                          Tuple{MOI.SingleVariable, MOI.Semicontinuous},
                                          Tuple{MOI.SingleVariable, MOI.Semiinteger},
                                          Tuple{MOI.VectorOfVariables, MOI.SOS1{T}},
                                          Tuple{MOI.VectorOfVariables, MOI.SOS2{T}}}

const LINEAR_FUNC_SET{T <: Real} = Union{Tuple{MOI.SingleVariable, MOI.LessThan{T}},
                                         Tuple{MOI.SingleVariable, MOI.GreaterThan{T}},
                                         Tuple{MOI.SingleVariable, MOI.EqualTo{T}},
                                         Tuple{MOI.SingleVariable, MOI.Interval{T}},
                                         Tuple{MOI.ScalarAffineFunction, MOI.LessThan{T}},
                                         Tuple{MOI.ScalarAffineFunction, MOI.GreaterThan{T}},
                                         Tuple{MOI.ScalarAffineFunction, MOI.EqualTo{T}},
                                         Tuple{MOI.ScalarAffineFunction, MOI.Interval{T}},
                                         Tuple{MOI.VectorOfVariables,    MOI.Nonnegatives},
                                         Tuple{MOI.VectorOfVariables,    MOI.Nonpositives},
                                         Tuple{MOI.VectorOfVariables,    MOI.Zeros},
                                         Tuple{MOI.VectorAffineFunction, MOI.Nonnegatives},
                                         Tuple{MOI.VectorAffineFunction, MOI.Nonpositives},
                                         Tuple{MOI.VectorAffineFunction, MOI.Zeros}}

const CONE_FUNC_SET =  Union{Tuple{MOI.VectorOfVariables, MOI.NormInfinityCone},
                             Tuple{MOI.VectorOfVariables, MOI.NormOneCone},
                             Tuple{MOI.VectorOfVariables, MOI.SecondOrderCone},
                             Tuple{MOI.VectorOfVariables, MOI.RotatedSecondOrderCone},
                             Tuple{MOI.VectorOfVariables, MOI.GeometricMeanCone},
                             Tuple{MOI.VectorOfVariables, MOI.ExponentialCone},
                             Tuple{MOI.VectorOfVariables, MOI.DualExponentialCone},
                             Tuple{MOI.VectorOfVariables, MOI.RelativeEntropyCone},
                             Tuple{MOI.VectorOfVariables, MOI.NormSpectralCone},
                             Tuple{MOI.VectorOfVariables, MOI.NormNuclearCone},
                             Tuple{MOI.VectorOfVariables, MOI.PowerCone},
                             Tuple{MOI.VectorOfVariables, MOI.DualPowerCone},
                             }

const SDP_FUNC_SET = Union{Tuple{MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeSquare},
                           Tuple{MOI.VectorOfVariables, MOI.RootDetConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.RootDetConeSquare},
                           Tuple{MOI.VectorOfVariables, MOI.LogDetConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.LogDetConeSquare}}

_in_prob(::Val{LP}, S, T)   = Tuple{S, T} <: LINEAR_FUNC_SET{Float64}
_in_prob(::Val{MILP}, S, T) = _in_prob(Val{LP}(), S, T) | (Tuple{S, T} <: INTEGER_FUNC_SET{Float64})
_in_prob(::Val{SOCP}, S, T) = _in_prob(Val{LP}(), S, T) | (Tuple{S, T} <: CONE_FUNC_SET)
_in_prob(::Val{SDP}, S, T)  = _in_prob(Val{SOCP}(), S, T) | (Tuple{S, T} <: SDP_FUNC_SET)

for F in (LP, MILP, SOCP, SDP)
    @eval function _parse_classify_problem(::Val{$F}, m::Optimizer)
        ip = _input_model(m)
        flag = isempty(filter(x -> !_in_prob(Val{$F}(), x[1], x[2]), MOI.get(ip, MOI.ListOfConstraints())))
        flag &= _input_nlp_data(m) === nothing
        if flag
            m._problem_type = $F
        end
        return flag
    end
end

#=
Labels branching and nonbranching variables, remakes the optimizer for correct
dimension of NodeBB.
=#
function _parse_classify_problem(::Val{MINCVX}, m::Optimizer{T}) where T<:AbstractFloat
    m._problem_type = MINCVX

    d = Dict{Int,Bool}()

    sqf_lt_list = MOI.get(m, MOI.ListOfConstraintIndices{SQF,LT}())
    sqf_gt_list = MOI.get(m, MOI.ListOfConstraintIndices{SQF,GT}())
    sqf_et_list = MOI.get(m, MOI.ListOfConstraintIndices{SQF,ET}())
    foreach(ci -> _label_branch!(d, m, ci), sqf_lt_list)
    foreach(ci -> _label_branch!(d, m, ci), sqf_gt_list)
    foreach(ci -> _label_branch!(d, m, ci), sqf_et_list)

    soc_list = MOI.get(m, MOI.ListOfConstraintIndices{VECVAR,SOC_CONE}())
    exp_list = MOI.get(m, MOI.ListOfConstraintIndices{VECVAR,EXP_CONE}())
    pow_list = MOI.get(m, MOI.ListOfConstraintIndices{VECVAR,POW_CONE}())
    psd_list = MOI.get(m, MOI.ListOfConstraintIndices{VECVAR,PSD_CONE}())
    foreach(ci -> _label_branch!(d, m, ci, start = 2), soc_list)
    foreach(ci -> _label_branch!(d, m, ci), exp_list)
    foreach(ci -> _label_branch!(d, m, ci), pow_list)
    foreach(ci -> _label_branch!(d, m, ci), psd_list)

    #Get branch variables from nonlinear expressions
    _label_nl!(d, _input_nlp_data(m))

    branch_num = length(keys(d))

    m._solver = GlobalOptimizer{branch_num, T, typeof(_ext_type(m))}()
    m._solver._ext_type = _ext_type(m)

    return true
end


"""
Classifies the problem type
"""
function _parse_classify_problem!(m::Optimizer)
    #_parse_bridge_nl!(m)
    _parse_classify_problem(Val{LP}(), m)       && return
    _parse_classify_problem(Val{MILP}(), m)     && return
    _parse_classify_problem(Val{SOCP}(), m)     && return
    _parse_classify_problem(Val{SDP}(), m)      && return
    _parse_classify_problem(Val{MINCVX}(), m)   && return
    return
end
