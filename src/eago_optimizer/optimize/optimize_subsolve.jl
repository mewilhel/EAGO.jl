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

for (T, optimizer_field) in ((LP, :lp_optimizer),
                             (MILP, :mip_optimizer),
                             (SOCP, :second_order_cone_optimizer),
                             (SDP, :semidefinite_optimizer),
                             (DIFF_CVX, :nlp_optimizer),
                             (MICVX, :minlp_optimizer))

    @eval function optimize!(::Val{$T}, m::Optimizer)
        opt = m.$optimizer_field
        MOI.empty!(opt)
        MOI.copy_to(opt, m._input_problem)
        set_config!(m, opt)

        if m._parameters.verbosity < 5
            MOI.set(opt, MOI.Silent(), true)
        end
        m._parse_time = time() - m._start_time

        MOI.optimize!(opt)
        _unpack_local_solve!(m, opt)
        return nothing
    end
end
