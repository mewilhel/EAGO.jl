# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eeago_optimizer/optimize/optimize_lp_cone.jl
# Contains the optimize! routines for LP, SOCP, (and in the future MILP and
# MISOCP) type problems. This also includes functions to add variables,
# linear constraints, soc constraints, and unpack solutions.
#############################################################################

#=
LP          -> COPY TO RELAXED SOLVER AND SOLVE
MILP        -> COPY TO RELAXED SOLVER AND SOLVE
SOCP        -> COPY TO RELAXED SOLVER AND SOLVE
MISOCP      -> COPY TO RELAXED SOLVER AND SOLVE
DIFF_CVX    -> COPY TO NLP SOLVER AND SOLVE (POTENTIAL MULTISTART)
DIFF_NCVX   -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
MINCVX      -> APPLY GLOBAL SOLVER (LOCAL SOLVE OPTION FUTURE FEATURE)
=#

function _unpack_final_subsolve!(m::Optimizer, opt::T; adjust_bnd::Bool = true) where T

    # TODO: ASSUMES NO BRIDGING NECESSARY, MAY NEED TO ADD SUPPORT LATER FOR BRIDGES
    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._primal_status_code = MOI.get(opt, MOI.PrimalStatus())

    if MOI.get(opt, MOI.ResultCount()) > 0
        variable_indices = [m._input_to_solution_map[vi] for vi in MOI.get(opt, MOI.ListOfVariableIndices())]
        m._solution = MOI.get(opt, MOI.VariablePrimal(), variable_indices)
        for (F, S) in MOI.get(m, MOI.ListOfConstraints())
            if !(F == SV)
                for ci in MOI.get(m, MOI.ListOfConstraintIndices{F,S}())
                    subsolver_ci = m._input_to_solution_map[ci]
                    primal_val = MOI.get(opt, MOI.ConstraintPrimal(), subsolver_ci)
                     _set_cons_primal!(m, ci, primal_val)
                 end
             end
        end
        m._objective_value = MOI.get(opt, MOI.ObjectiveValue())
        m._objective_bound = MOI.get(opt, MOI.ObjectiveBound())
    end

    return nothing
end

for S in (LP, MILP, SOCP, SDP)
    @eval _bridge_optimizer!(v::Val{$S}, m) = nothing
end
for S in (DIFF_CVX, MICVX)
    @eval function _bridge_optimizer!(v::Val{$S}, m)
        MOIB.add_bridge(m, MOIB.Constraint.SOCtoNonConvexQuadBridge{Float64})
        MOIB.add_bridge(m, MOIB.Constraint.RSOCtoNonConvexQuadBridge{Float64})
        return
    end
end

for (T, sub_optimizer) in ((LP, :_lp_optimizer),
                            (MILP, :_mip_optimizer),
                            (SOCP, :_socp_optimizer),
                            (SDP, :_semidefinite_optimizer),
                            (DIFF_CVX, :_nlp_optimizer),
                            (MICVX, :_minlp_optimizer))

    @eval function optimize!(::Val{$T}, m::Optimizer)

        opt = MOI.instantiate(($sub_optimizer)(m), with_bridge_type = Float64)
        #set_config!(m, opt)
        _bridge_optimizer!(Val{$T}(), opt)
        m._input_to_solution_map = MOI.copy_to(opt,
                                               _input_model(m);
                                               copy_names = false)

        if _verbosity(m) < 5
            MOI.set(opt, MOI.Silent(), true)
        end
        m._parse_time = time() - m._start_time

        MOI.optimize!(opt)
        _unpack_final_subsolve!(m, opt)
    end
end
