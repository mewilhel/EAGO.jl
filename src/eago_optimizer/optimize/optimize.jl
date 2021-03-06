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
include("nonconvex/presolve.jl")
include("nonconvex/branch/branch.jl")
include("nonconvex/optimize_local_nlp.jl")
include("nonconvex/optimize_relax_lp.jl")
include("nonconvex/preprocess.jl")
include("nonconvex/postprocess.jl")
include("nonconvex/optimize_nonconvex.jl")

throw_optimize_hook!(m::Optimizer) = optimize_hook!(_ext_type(m), m)

function MOI.optimize!(m::Optimizer)

    m._start_time = time()
    m._time_left = m._time_limit

    if !_enable_optimize_hook(m)

        _parse_classify_problem!(m)

        new_time = time() - m._start_time
        m._parse_time = new_time
        m._run_time = new_time

        optimize!(Val{m._problem_type}(), m)
    else
        throw_optimize_hook!(m)
    end

    return nothing
end
