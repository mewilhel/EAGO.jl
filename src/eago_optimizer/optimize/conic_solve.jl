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
NS_CVX      -> COPY TO NLP SOLVER AND SOLVE (POTENTIAL MULTISTART)
DIFF_NCVX   -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
NS_NCVX     -> APPLY GLOBAL SOLVER (UNLESS USER REQUEST LOCAL SOLVE THEN NLP)
MINCVX      -> APPLY GLOBAL SOLVER (LOCAL SOLVE OPTION FUTURE FEATURE)
=#

### LP and MILP routines
function _add_soc_constraints!(m::Optimizer, opt::T) where T
    if MOI.supports_constraint(opt, VECOFVAR, SECOND_ORDER_CONE)
        ip = m._input_problem
        foreach(fs -> MOI.add_constraint!(opt, fs[1], fs[2]), _second_order_cone_constraints(ip))
    else
        # TODO: FIX
        for (func, set) in m._input_problem._conic_second_order
            quad_terms = SQT[SQT((), func.variables[i], func.variables[i]) for i = 1:length(func.variables)]
            sqf = SQF(SQT[], SAF[], 0.0)
            MOI.add_constraint(opt, sqf, LT_ZERO)
        end
    end
    return nothing
end

function optimize!(::Val{SOCP}, m::Optimizer)

    socp_optimizer = m.second_order_cone_optimizer
    set_default_config!(socp_optimizer)
    MOI.empty!(socp_optimizer)

    ip = m._input_problem
    m._relaxed_variable_index = add_variables(m, socp_optimizer, _variable_count(ip))
    _add_linear_constraints!(ip, socp_optimizer)
    _add_soc_constraints!(ip, socp_optimizer)
    _add_sv_or_aff_obj!(ip, socp_optimizer)
    _set_sense!(ip, socp_optimizer)

    if m._parameters.verbosity < 5
        MOI.set(socp_optimizer, MOI.Silent(), true)
    end
    m._parse_time = time() - m._start_time

    MOI.optimize!(socp_optimizer)
    _unpack_local_solve!(m, socp_optimizer)

    return nothing
end

optimize!(::Val{MISOCP}, m::Optimizer) = optimize!(Val{SOCP}(), m)
