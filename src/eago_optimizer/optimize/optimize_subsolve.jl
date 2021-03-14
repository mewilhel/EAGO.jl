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

function _unpack_final_solve!(m::Optimizer, opt::T, idx_map; adjust_bnd::Bool = true) where T

    m._termination_status_code = MOI.get(opt, MOI.TerminationStatus())
    m._result_status_code = MOI.get(opt, MOI.PrimalStatus())

    if MOI.get(opt, MOI.ResultCount()) > 0
        variable_indices = MOI.get(opt, MOI.ListOfVariableIndices())
        m._solution = MOI.get(opt, MOI.VariablePrimal(), variable_indices)
        for (F, S) = MOI.get(opt, MOI.ListOfConstraints())
            for ci in MOI.ListOfConstraintIndices{F,S}()
                primal_val = MOI.get(opt, MOI.ConstraintPrimal(), ci)
                 _set_cons_primal!(m, ci, primal_val)
             end
        end
        m._objective_value = MOI.get(opt, MOI.ObjectiveValue())
        m._objective_bound = MOI.get(opt, MOI.ObjectiveBound())
    end

    return nothing
end

for (T, optimizer_field) in ((LP, :lp_optimizer),
                             (MILP, :mip_optimizer),
                             (SOCP, :second_order_cone_optimizer),
                             (SDP, :semidefinite_optimizer),
                             (DIFF_CVX, :nlp_optimizer),
                             (MICVX, :minlp_optimizer))

    @eval function optimize!(::Val{$T}, m::Optimizer)

        opt = MOIB.full_bridge_optimizer(m.$optimizer_field, T)
        MOI.empty!(opt)
        idx_map = MOI.copy_to(opt, m._input_problem)
        #set_config!(m, opt)

        if m.verbosity < 5
            MOI.set(opt, MOI.Silent(), true)
        end
        m._parse_time = time() - m._start_time

        MOI.optimize!(opt)
        _unpack_final_solve!(m, opt, idx_map)
    end
end
