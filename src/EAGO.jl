# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/EAGO.jl
# The main file for EAGO.
#############################################################################

__precompile__()

module EAGO

    import MathOptInterface

    using Reexport, Cassette, IntervalArithmetic, NumericIO, DocStringExtensions
    using FastRounding, SpecialFunctions

    import JuMP
    import JuMP._Derivatives: operators, NodeData, VARIABLE, VALUE, SUBEXPRESSION,
                              Linearity, classify_linearity, UserOperatorRegistry,
                              register_multivariate_operator!
    using JuMP._Derivatives: univariate_operators,
                             univariate_operator_to_id

    import Ipopt, Hypatia, GLPK, Pavito

    using DataStructures: OrderedDict, BinaryMinMaxHeap, popmin!, popmax!, top
    using SparseArrays: SparseMatrixCSC, spzeros, rowvals, nzrange, nonzeros, sparse, findnz
    using LinearAlgebra: eigmin, norm

    import IntervalArithmetic: mid

    @reexport using McCormick
    @reexport using SpecialFunctions

    using IntervalContractors

    const MOI = MathOptInterface
    const MOIU = MathOptInterface.Utilities
    const MOIB = MathOptInterface.Bridges

    const SV = MOI.SingleVariable
    const SAF = MOI.ScalarAffineFunction{Float64}
    const SQF = MOI.ScalarQuadraticFunction{Float64}
    const VECVAR = MOI.VectorOfVariables

    const SQT = MOI.ScalarQuadraticTerm{Float64}
    const SAT = MOI.ScalarAffineTerm{Float64}

    const LT = MOI.LessThan{Float64}
    const GT = MOI.GreaterThan{Float64}
    const ET = MOI.EqualTo{Float64}
    const IT = MOI.Interval{Float64}
    const ZO = MOI.ZeroOne
    const SOC_CONE = MOI.SecondOrderCone
    const POW_CONE = MOI.PowerCone
    const EXP_CONE = MOI.ExponentialCone
    const PSD_CONE = MOI.PositiveSemidefiniteConeTriangle

    const VI = MOI.VariableIndex
    const CI = MOI.ConstraintIndex
    const ListCI = MOI.ListOfConstraintIndices

    const SCoefC = MOI.ScalarCoefficientChange
    const SConsC = MOI.ScalarConstantChange
    const MOIU = MOI.Utilities

    const LT_ZERO = LT(0.0)

    # Add storage types for EAGO optimizers
    export NodeBB, get_history, get_lower_bound, get_upper_bound, get_lower_time,
           get_upper_time, get_preprocess_time, get_postprocess_time, get_lower_bound, get_solution_time,
           get_iteration_number, get_node_count, get_absolute_gap, get_relative_gap

    export register_eago_operators!

    # map/reduce non-allocating no bounds checking map-reduce like utilities
    include("eago_optimizer/utilities/unsafe_utilities.jl")
    include("eago_optimizer/utilities/register_special.jl")

    # defines structure used to store information at each iteration of global optimize
    include("eago_optimizer/logging/log.jl")

    # defines the optimizer structures
    include("eago_optimizer/types/types.jl")

    # defines routines to add saf, sqf, and nlp block constraints
    include("eago_optimizer/math_opt_interface.jl")

    # functions which print information to console
    include("eago_optimizer/display.jl")

    #
    include("eago_optimizer/evaluate/evaluate.jl")

    #
    include("eago_optimizer/domain_reduction.jl")

    #
    include("eago_optimizer/parse_classify.jl")

    #
    include("eago_optimizer/logging/log_iteration.jl")

    #
    include("eago_optimizer/optimize/optimize.jl")

    # import the script solving utilities
    include("eago_script/script.jl")

    # routines for solving SIPs
    export SIPResult, SIPProblem, SIPCallback, SIPSubResult,
           sip_solve, SIPRes, SIPResRev, SIPHybrid,
           build_model, set_tolerance_inner!, set_tolerance!, get_disc_set,
           sip_llp!, sip_bnd!, sip_res!, get_sip_optimizer, check_convergence,
           LowerLevel1, LowerLevel2, LowerLevel3, LowerProblem, UpperProblem,
           ResProblem, AbstractSIPAlgo, AbstractSubproblemType
    include("eago_semiinfinite/semiinfinite.jl")
end
