# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/optimize.jl
# Load subproblem type optimization routines and define general optimize!
# throw_optimize_hook!.
#############################################################################

include("optimize_subsolve.jl")
include("optimize_convex.jl")
include("optimize_nonconvex.jl")

throw_optimize_hook!(m::Optimizer) = optimize_hook!(m.ext_type, m)

function _setup_timers!(m::Optimizer)
    m._start_time = time()
    m._time_left = m.time_limit
    return
end
function _set_parse_time!(m::Optimizer)
    new_time = time() - m._start_time
    m._parse_time = new_time
    m._run_time = new_time
    return
end


function MOI.optimize!(m::Optimizer)

    _setup_timers!(m)

    if !m.enable_optimize_hook
        _parse_classify_problem!(m)
        _set_parse_time!(m)
        optimize!(Val{m._working_problem._problem_type}(), m)
    else
        throw_optimize_hook!(m)
    end

    return nothing
end
