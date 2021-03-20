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
                             Tuple{MOI.VectorOfVariables, MOI.NormNuclearCone}}

const SDP_FUNC_SET = Union{Tuple{MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeSquare},
                           Tuple{MOI.VectorOfVariables, MOI.RootDetConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.RootDetConeSquare},
                           Tuple{MOI.VectorOfVariables, MOI.LogDetConeTriangle},
                           Tuple{MOI.VectorOfVariables, MOI.LogDetConeSquare}}

_in_prob(::Val{LP}, S, T)   = Tuple{S, T} <: LINEAR_FUNC_SET{Float64}
_in_prob(::Val{MILP}, S, T) = _in_linear(S, T) | (Tuple{S, T} <: INTEGER_FUNC_SET{Float64})
_in_prob(::Val{SOCP}, S, T) = _in_linear(S, T) | (Tuple{S, T} <: CONE_FUNC_SET)
_in_prob(::Val{SDP}, S, T)  = _in_socp(S, T) | (Tuple{S, T} <: SDP_FUNC_SET)

for F in (LP, MILP, SOCP, SDP)
    @eval function _parse_classify_problem(::Val{LP}, m::Optimizer)
        ip = _input_model(m)
        flag &= isempty(filter(x -> !_in_prob(Val{$F}(), x[1], x[2]), MOI.get(ip, MOI.ListOfConstraints())))
        flag &= _input_nlp_data(m) === nothing
        if flag
            m._problem_type = $F
        end
        return flag
    end
end

function _parse_classify_problem(::Val{MINCVX}, m::Optimizer)
    m._problem_type = MINCVX
    return true
end


"""
Classifies the problem type
"""
function _parse_classify_problem!(m::Optimizer)
    _parse_classify_problem(Val{LP}(), m)       && return
    _parse_classify_problem(Val{MILP}(), m)     && return
    _parse_classify_problem(Val{SOCP}(), m)     && return
    _parse_classify_problem(Val{SDP}(), m)      && return
    _parse_classify_problem(Val{MINCVX}(), m)   && return
    return
end
