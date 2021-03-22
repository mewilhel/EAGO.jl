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
function node_selection!(t::ExtensionType, m::GlobalOptimizer{N,T}) where {N,T}
    m._node_count -= 1
    m._current_node = popmin!(m._stack)
    return nothing
end

"""
$(SIGNATURES)

Stores the current node to the stack after updating lower/upper bounds.
"""
function single_storage!(t::ExtensionType, m::GlobalOptimizer{N,T})
    y = m._current_node
    m._node_repetitions += 1
    m._node_count += 1
    lower_bound = max(y.lower_bound, m._lower_objective_value)
    upper_bound = min(y.upper_bound, m._upper_objective_value)
    push!(m._stack, NodeBB(y.lower_variable_bound, y.upper_variable_bound,
                           lower_bound, upper_bound, y.depth, y.id))

    return nothing
end

"""
$(SIGNATURES)

Selects and deletes nodes from stack with lower bounds greater than global
upper bound.
"""
function fathom!(t::ExtensionType, m::GlobalOptimizer{N,T})

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

"""
$(SIGNATURES)

Checks for termination of algorithm due to satisfying absolute or relative
tolerance, infeasibility, or a specified limit, returns a boolean valued true
if algorithm should continue.
"""
function termination_check(t::ExtensionType, m::GlobalOptimizer{N,T}) where {N,T<:Real}

    node_in_stack = length(m._stack)
    L = m._global_lower_bound
    U = m._global_upper_bound

    if node_in_stack == 0

        if m._first_solution_node > 0
            m._termination_status_code = MOI.OPTIMAL
            m._result_status_code = MOI.FEASIBLE_POINT
            (m.verbosity >= 3) && println("Empty Stack: Exhaustive Search Finished")

        else
            m._termination_status_code = MOI.INFEASIBLE
            m._result_status_code = MOI.INFEASIBILITY_CERTIFICATE
            (m.verbosity >= 3) && println("Empty Stack: Infeasible")
        end

    elseif node_in_stack >= m.node_limit

        m._termination_status_code = MOI.NODE_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m.verbosity >= 3) && println("Node Limit Exceeded")

    elseif m._iteration_count >= m.iteration_limit

        m._termination_status_code = MOI.ITERATION_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m.verbosity >= 3) && println("Maximum Iteration Exceeded")

    elseif ~relative_tolerance(L, U, m.relative_tolerance)

        m._termination_status_code = MOI.OPTIMAL
        m._result_status_code = MOI.FEASIBLE_POINT
        (m.verbosity >= 3) && println("Relative Tolerance Achieved")

    elseif (U - L) < m.absolute_tolerance

        m._termination_status_code = MOI.OPTIMAL
        m._result_status_code = MOI.FEASIBLE_POINT
        (m.verbosity >= 3) && println("Absolute Tolerance Achieved")

    elseif m._run_time > m.time_limit

        m._termination_status_code = MOI.TIME_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m.verbosity >= 3) && println("Time Limit Exceeded")

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
function convergence_check(t::ExtensionType, m::GlobalOptimizer{N,T}) where {N,T<:Real}

  L = m._lower_objective_value
  U = m._global_upper_bound
  t = (U - L) <= _absolute_tolerance(m)
  if (U < Inf) && (L > Inf)
      t |= (abs(U - L)/(max(abs(L), abs(U))) <= _relative_tolerance(m)
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
optimize_hook!(t::ExtensionType, m::GlobalOptimizer{N,T}) where {N,T<:Real} = nothing

function store_candidate_solution!(m::GlobalOptimizer{N,T}) where {N,T<:Real}
    if m._upper_feasibility && (m._upper_objective_value < m._global_upper_bound)
        m._feasible_solution_found = true
        m._first_solution_node = m._maximum_node_id
        m._global_upper_bound = m._upper_objective_value
        @__dot__ m._solution = m._upper_solution
    end
    return nothing
end

function set_global_lower_bound!(m::GlobalOptimizer{N,T}) where {N,T<:Real}
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
parse_global!(m::GlobalOptimizer{N,T}) where T = parse_global!(m.ext_type, m)
presolve_global!(m::GlobalOptimizer{N,T}) where T = presolve_global!(m.ext_type, m)
termination_check(m::GlobalOptimizer{N,T}) where T = termination_check(m.ext_type, m)
cut_condition(m::GlobalOptimizer{N,T}) where T = cut_condition(m.ext_type, m)
convergence_check(m::GlobalOptimizer{N,T}) where T = convergence_check(m.ext_type, m)
repeat_check(m::GlobalOptimizer{N,T}) where T = repeat_check(m.ext_type, m)
node_selection!(m::GlobalOptimizer{N,T}) where T = node_selection!(m.ext_type, m)
preprocess!(m::GlobalOptimizer{N,T}) where T = preprocess!(m.ext_type, m)
lower_problem!(m::GlobalOptimizer{N,T}) where T = lower_problem!(m.ext_type, m)
add_cut!(m::GlobalOptimizer{N,T}) where T = add_cut!(m.ext_type, m)
upper_problem!(m::GlobalOptimizer{N,T}) where T = upper_problem!(m.ext_type, m)
postprocess!(m::GlobalOptimizer{N,T}) where T = postprocess!(m.ext_type, m)
single_storage!(m::GlobalOptimizer{N,T}) where T = single_storage!(m.ext_type, m)
branch_node!(m::GlobalOptimizer{N,T}) where T = branch_node!(m.ext_type, m)
fathom!(m::GlobalOptimizer{N,T}) where T = fathom!(m.ext_type, m)
revert_adjusted_upper_bound!(m::GlobalOptimizer{N,T}) where T = revert_adjusted_upper_bound!(m.ext_type, m)

"""
$(TYPEDSIGNATURES)

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(m::GlobalOptimizer)

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
