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

function _unpack_final_subsolve!(m::Optimizer, opt::T, idx_map; adjust_bnd::Bool = true) where T

    # TODO: ASSUMES NO BRIDGING NECESSARY, MAY NEED TO ADD SUPPORT LATER FOR BRIDGES
    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())

    if MOI.get(opt, MOI.ResultCount()) > 0
        variable_indices = [idx_map[vi] for vi in MOI.get(opt, MOI.ListOfVariableIndices())]
        m._solution = MOI.get(opt, MOI.VariablePrimal(), variable_indices)
        for (F, S) in MOI.get(m, MOI.ListOfConstraints())
            if !(F == SV)
                for ci in MOI.get(m, MOI.ListOfConstraintIndices{F,S}())
                    subsolver_ci = idx_map[ci]
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

_bridge_optimizer(::Val{LP}, m) = MOIB.full_bridge_optimizer(m, Float64)
_bridge_optimizer(::Val{MILP}, m) = MOIB.full_bridge_optimizer(m, Float64)
_bridge_optimizer(::Val{SOCP}, m) = MOIB.full_bridge_optimizer(m, Float64)
_bridge_optimizer(::Val{SDP}, m) = MOIB.full_bridge_optimizer(m, Float64)

for S in (DIFF_CVX, MICVX)
    @eval function _bridge_optimizer(::Val{$S}, m)
        opt = MOIB.full_bridge_optimizer(m, Float64)
        MOIB.add_bridge(opt, MOIB.Constraint.SOCtoNonConvexQuadBridge{Float64})
        MOIB.add_bridge(opt, MOIB.Constraint.RSOCtoNonConvexQuadBridge{Float64})
        return opt
    end
end

for (T, optimizer_field) in ((LP, :lp_optimizer),
                             (MILP, :mip_optimizer),
                             (SOCP, :socp_optimizer),
                             (SDP, :semidefinite_optimizer),
                             (DIFF_CVX, :nlp_optimizer),
                             (MICVX, :minlp_optimizer))

    @eval function optimize!(::Val{$T}, m::Optimizer)

        opt = m.$optimizer_field
        #set_config!(m, opt)
        bridged_opt = _bridge_optimizer(Val{$T}(), opt)
        idx_map = MOIU.automatic_copy_to(bridged_opt, m._input_problem; copy_names = false)

        if m.verbosity < 5
            MOI.set(bridged_opt, MOI.Silent(), true)
        end
        m._parse_time = time() - m._start_time

        MOI.optimize!(bridged_opt)
        _unpack_final_subsolve!(m, bridged_opt, idx_map)
    end
end
