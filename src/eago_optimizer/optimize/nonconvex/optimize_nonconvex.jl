# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/optimize/optimize_nonconvex.jl
# Contains the optimize! routine and subroutines needed in the branch and
# bound routine called by EAGO.
#############################################################################

"""
$(SIGNATURES)

Selects node with the lowest lower bound in stack.
"""
function node_selection!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    m._node_count -= 1
    m._current_node = popmin!(m._stack)
    return nothing
end

"""
$(SIGNATURES)

Stores the current node to the stack after updating lower/upper bounds.
"""
function single_storage!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    n = m._current_node
    m._node_repetitions += 1
    m._node_count += 1
    lower_bound = max(n.lower_bound, m._lower_objective_value)
    upper_bound = min(n.upper_bound, m._upper_objective_value)
    push!(m._stack, NodeBB(n.lower_variable_bound, n.upper_variable_bound,
                           lower_bound, upper_bound, n.depth, n.id))

    return nothing
end

"""
$(SIGNATURES)

Selects and deletes nodes from stack with lower bounds greater than global
upper bound.
"""
function fathom!(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    upper = m._global_upper_bound
    continue_flag = !isempty(m._stack)
    while continue_flag
        max_node = maximum(m._stack)
        max_check = (max_node.lower_bound > upper)
        if max_check
            popmax!(m._stack)
            m._node_count -= 1
            if isempty(m._stack)
                continue_flag = false
            end
        else
            if !max_check
                continue_flag = false
            elseif isempty(m._stack)
                continue_flag = false
            end
        end
    end
    return nothing
end

"""
$(SIGNATURES)

Checks to see if current node should be reprocessed.
"""
repeat_check(t::ExtensionType, m::GlobalOptimizer) = false

relative_gap(L::Float64, U::Float64) = ((L > -Inf) && (U < Inf)) ?  abs(U - L)/(max(abs(L), abs(U))) : Inf
relative_tolerance(L::Float64, U::Float64, tol::Float64) = relative_gap(L, U)  > tol || ~(L > -Inf)

function _terminate_store!(m::GlobalOptimizer, t::MOI.TerminationStatusCode,
                           r::MOI.ResultStatusCode, s::String)
    m._termination_status_code = t
    m._result_status_code = r
    (_verbosity(m) >= 3) && println(s)
    return nothing
end

"""
$(SIGNATURES)

Checks for termination of algorithm due to satisfying absolute or relative
tolerance, infeasibility, or a specified limit, returns a boolean valued true
if algorithm should continue.
"""
function termination_check(t::ExtensionType, m::GlobalOptimizer)
    if length(m._stack) == 0
        if m._first_solution_node > 0
            _terminate_store!(m, MOI.OPTIMAL, MOI.FEASIBLE_POINT, "Empty Stack: Exhaustive Search Finished")
        else
            _terminate_store!(m, MOI.INFEASIBLE,MOI.INFEASIBILITY_CERTIFICATE, "Empty Stack: Infeasible")
        end
    elseif length(m._stack) >= m._node_limit
        _terminate_store!(m, MOI.NODE_LIMIT, MOI.UNKNOWN_RESULT_STATUS, "Node Limit Exceeded")
    elseif m._iteration_count >= m.iteration_limit
        _terminate_store!(m, MOI.ITERATION_LIMIT, MOI.UNKNOWN_RESULT_STATUS, "Maximum Iteration Exceeded")
    elseif _relative_gap(m) < m._relative_tolerance
        _terminate_store!(m, MOI.OPTIMAL, MOI.FEASIBLE_POINT, "Relative Tolerance Achieved")
    elseif _relative_gap(m) < m._absolute_tolerance
        _terminate_store!(m, MOI.OPTIMAL, MOI.FEASIBLE_POINT, "Absolute Tolerance Achieved")
    elseif m._run_time > m.time_limit
        _terminate_store!(m, MOI.TIME_LIMIT, MOI.UNKNOWN_RESULT_STATUS, "Time Limit Exceeded")
    else
        return false
    end
    return true
end

"""
$(SIGNATURES)

Checks for convergence of algorithm with respect to absolute and/or relative
tolerances.
"""
function convergence_check(t::ExtensionType, m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

  L = m._lower_objective_value
  U = m._global_upper_bound
  t = (U - L) <= _absolute_tolerance(m)
  if (U < Inf) && (L > Inf)
      t |= (abs(U - L)/(max(abs(L), abs(U))) <= _relative_tolerance(m))
  end

  if t && m._min_converged_value < Inf
      m._min_converged_value = min(m._min_converged_value, L)
  else
      m._min_converged_value = L
  end

  return t
end

"""
$(SIGNATURES)

Provides a hook for extensions to EAGO as opposed to standard global, local,
or linear solvers.
"""
optimize_hook!(t::ExtensionType, m::GlobalOptimizer) = nothing

function store_candidate_solution!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    if m._upper_feasibility && (m._upper_objective_value < m._global_upper_bound)
        m._feasible_solution_found = true
        m._first_solution_node = m._maximum_node_id
        m._global_upper_bound = m._upper_objective_value
        @__dot__ m._solution = m._upper_solution
    end
    return nothing
end

function set_global_lower_bound!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
    if !isempty(m._stack)
        min_node = minimum(m._stack)
        lower_bound = min_node.lower_bound
        if m._global_lower_bound < lower_bound
            m._global_lower_bound = lower_bound
        end
    end
    return nothing
end

# wraps subroutine call to isolate ExtensionType
for f in (:parse_global!, :presolve_global!, :termination_check, :cut_condition,
          :convergence_check, :repeat_check, :node_selection!, :preprocess!,
          :lower_problem!, :add_cut!, :upper_problem!, :postprocess!,
          :single_storage, :branch_node!, :fathom!)
    @eval function ($f)(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}
        ($f)(m.ext_type, m)
    end
end

"""
$(TYPEDSIGNATURES)

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(m::GlobalOptimizer{N,T,S}) where {N,T<:AbstractFloat,S}

    m._iteration_count = 1
    m._node_count = 1

    parse_global!(m)
    presolve_global!(m)

    logging_on = m.log_on
    verbosity = m.verbosity

    # terminates when max nodes or iteration is reach, or when node stack is empty
    while !termination_check(m)

        # Selects node, deletes it from stack, prints based on verbosity
        node_selection!(m)
        (verbosity >= 3) && print_node!(m)

        # Performs prepocessing and times
        logging_on && (start_time = time())
        preprocess!(m)
        if logging_on
            m._last_preprocess_time = time() - start_time
        end

        if m._preprocess_feasibility

            # solves & times lower bounding problem
            logging_on && (start_time = time())
            m._cut_iterations = 1
            lower_problem!(m)
            while cut_condition(m)
                add_cut!(m)
            end
            if logging_on
                m._last_lower_problem_time = time() - start_time
            end
            print_results!(m, true)
            print_results_post_cut!(m)

            # checks for infeasibility stores solution
            if m._lower_feasibility
                if !convergence_check(m)

                    logging_on && (start_time = time())
                    upper_problem!(m)
                    if logging_on
                        m._last_upper_problem_time = time() - start_time
                    end
                    print_results!(m, false)
                    store_candidate_solution!(m)
                    if m._input_problem._optimization_sense === MOI.FEASIBILITY_SENSE
                        if !m.feasible_local_continue || m.local_solve_only
                            break
                        end
                    end

                    # Performs and times post processing
                    logging_on && (start_time = time())
                    postprocess!(m)
                    if logging_on
                        m._last_postprocessing_time = time() - start_time
                    end

                    # Checks to see if the node
                    if m._postprocess_feasibility
                        if repeat_check(m)
                            single_storage!(m)
                        else
                            branch_node!(m)
                        end
                    end
                end
            end
            fathom!(m)
        else
            m._lower_objective_value = -Inf
            m._lower_feasibility = false
            m._upper_feasibility = false
        end
        set_global_lower_bound!(m)
        m._run_time = time() - m._start_time
        m._time_left = m.time_limit - m._run_time
        log_iteration!(m)
        print_iteration!(m)
        m._iteration_count += 1
    end

    m._objective_value = m._global_upper_bound

    # Prints the solution
    print_solution!(m)

    return nothing
end

optimize!(::Val{MINCVX}, m::Optimizer) = global_solve!(m._solver)
