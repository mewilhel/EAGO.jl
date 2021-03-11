# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/evaluator/passes.jl
# Functions used to compute forward pass of nonlinear functions which include:
# set_value_post, overwrite_or_intersect, forward_pass_kernel, associated blocks
#############################################################################

const FORWARD_DEBUG = false

"""
    set_value_post

Post process set_value operator. By default, performs the affine interval cut on
a MC structure.
"""
function set_value_post(x_values::Vector{Float64}, val::MC{N,T}, lower_variable_bounds::Vector{Float64},
                        upper_variable_bounds::Vector{Float64}, sparsity::Vector{Int},
                        subgrad_tol::Float64) where {N, T<:RelaxTag}

    lower = val.cv
    upper = val.cc
    lower_refinement = true
    upper_refinement = true

    for i = 1:N

        cv_val = @inbounds val.cv_grad[i]
        cc_val = @inbounds val.cc_grad[i]

        i_sol = @inbounds sparsity[i]
        x_val = @inbounds x_values[i_sol]
        lower_bound = @inbounds lower_variable_bounds[i_sol]
        upper_bound = @inbounds upper_variable_bounds[i_sol]

        if lower_refinement
            if cv_val > 0.0
                if isinf(lower_bound)
                    !upper_refinement && break
                    lower_refinement = false
                else
                    lower += cv_val*(lower_bound - x_val)
                    #delX = sub_round(lower_bound, x_val, RoundDown)
                    #lower = add_round(lower, mul_round(cv_val, delX, RoundDown), RoundDown)
                end
            else
                if isinf(upper_bound)
                    !upper_refinement && break
                    lower_refinement = false
                else
                    lower += cv_val*(upper_bound - x_val)
                    #delX = sub_round(upper_bound, x_val, RoundUp)
                    #lower = add_round(lower, mul_round(cv_val, delX, RoundDown), RoundDown)
                end
            end
        end

        if upper_refinement
            if cc_val > 0.0
                if isinf(lower_bound)
                    !lower_refinement && break
                    upper_refinement = false
                else
                    upper += cc_val*(upper_bound - x_val)
                    #delX = sub_round(upper_bound, x_val, RoundUp)
                    #upper = add_round(upper, mul_round(cc_val, delX, RoundUp), RoundUp)
                end
            else
                if isinf(upper_bound)
                    !lower_refinement && break
                    upper_refinement = false
                else
                    upper += cc_val*(lower_bound - x_val)
                    #delX = sub_round(lower_bound, x_val, RoundDown)
                    #upper = add_round(upper, mul_round(cc_val, delX, RoundUp), RoundUp)
                end
            end
        end
    end

    if lower_refinement && (val.Intv.lo + subgrad_tol > lower)
        lower = val.Intv.lo
    elseif !lower_refinement
        lower = val.Intv.lo
    else
        lower -= subgrad_tol #sub_round(lower, subgrad_tol, RoundDown)
    end

    if upper_refinement && (val.Intv.hi - subgrad_tol < upper)
        upper = val.Intv.hi
    elseif !upper_refinement
        upper = val.Intv.hi
    else
        upper += subgrad_tol #add_round(upper, subgrad_tol, RoundUp)
    end

    return MC{N,T}(val.cv, val.cc, Interval{Float64}(lower, upper), val.cv_grad, val.cc_grad, val.cnst)
end

"""
    overwrite_or_intersect

Intersects the new set valued operator with the prior and performs affine bound tightening

- First forward pass: `is_post` should be set by user option, `is_intersect` should be false
  so that the tape overwrites existing values, and the `interval_intersect` flag could be set
  to either value.
- Forward CP pass (assumes same reference point): `is_post` should be set by user option,
  `is_intersect` should be true so that the tape intersects with  existing values, and the
  `interval_intersect` flag should be false.
- Forward CP pass (assumes same reference point): `is_post` should be set by user option,
  `is_intersect` should be true so that the tape intersects with existing values, and the
  `interval_intersect` flag should be false.
- Subsequent forward passes at new points: is_post` should be set by user option,
  `is_intersect` should be true so that the tape intersects with existing values, and the
  `interval_intersect` flag should be `true` as predetermined interval bounds are valid but
   the prior values may correspond to different points of evaluation.
"""
function overwrite_or_intersect(xMC::MC{N,T}, past_xMC::MC{N,T}, x::Vector{Float64}, lbd::Vector{Float64},
                                ubd::Vector{Float64}, subgrad_tol::Float64, sparsity::Vector{Int}, is_post::Bool,
                                is_intersect::Bool,
                                interval_intersect::Bool) where {N,T<:RelaxTag}

    if is_post && is_intersect && interval_intersect
        return set_value_post(x, xMC ∩ past_xMC.Intv, lbd, ubd, sparsity, subgrad_tol)

    elseif is_post && is_intersect && !interval_intersect
        return set_value_post(x, xMC ∩ past_xMC, lbd, ubd, sparsity, subgrad_tol)

    elseif is_post && !is_intersect
        return set_value_post(x, xMC, lbd, ubd, sparsity, subgrad_tol)

    elseif !is_post && is_intersect && interval_intersect
        return xMC ∩ past_xMC.Intv

    elseif !is_post && is_intersect && !interval_intersect
        return xMC ∩ past_xMC

    end
    return xMC
end

"""
    forward_plus_binary!

Updates storage tapes with forward evalution of node representing `n = x + y`.
"""
function forward_plus_binary!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                              numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                              x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, subgrad_tol::Float64,
                              sparsity::Vector{Int},
                              is_post::Bool, is_intersect::Bool, is_first_eval::Bool, interval_intersect::Bool) where {N,T<:RelaxTag}

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index =  children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = setstorage[arg2_index]
        num2 = 0.0
    end

    output_is_number = arg1_is_number && arg2_is_number

    # a + b
    if output_is_number
         numberstorage[k] = num1 + num2

    # x + b
    elseif !arg1_is_number && arg2_is_number
        outset = set1 + num2
        # is_first_eval ? (set1 + num2) : plus_kernel(set1, num2, setstorage[k].Intv)

    # a + y
    elseif arg1_is_number && !arg2_is_number
        outset = num1 + set2
        # is_first_eval ? (num1 + set2) : plus_kernel(num1, set2, setstorage[k].Intv)

    # x + y
    else
        outset = set1 + set2
        # is_first_eval ? (set1 + set2) : plus_kernel(set1, set2, setstorage[k].Intv)

    end

     numvalued[k] = output_is_number
    if !output_is_number
         setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, false,
                                                         is_intersect, interval_intersect)
    end

    return nothing
end

"""
    forward_plus_narity!

Updates storage tapes with forward evalution of node representing `n = +(x, y, z,...)`.
"""
function forward_plus_narity!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                              numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                              x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64},
                              subgrad_tol::Float64, sparsity::Vector{Int},
                              is_post::Bool, is_intersect::Bool, interval_intersect::Bool) where {N,T<:RelaxTag}


    # get row indices
    idx = first(children_idx)

    # extract values for argument 1
    arg_index =  children_arr[idx]
    output_is_number = numvalued[arg_index]
    if output_is_number
       tmp_set = zero(MC{N,T})
       tmp_num = numberstorage[arg_index]
    else
       tmp_num = 0.0
       tmp_set = setstorage[arg_index]
    end
    output_is_number = true

    for idx = 2:length(children_idx)
        cidx = children_idx[idx]
        arg_index =  children_arr[cidx]
        arg_is_number = numvalued[arg_index]
        if arg_is_number
            tmp_num += numberstorage[arg_index]
        else
            tmp_set += setstorage[arg_index]
        end
        output_is_number &= arg_is_number
    end

     numvalued[k] = output_is_number
    if output_is_number
         numberstorage[k] = tmp_num
    else
        tmp_set += tmp_num
        setstorage[k] = overwrite_or_intersect(tmp_set, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, false, is_intersect,
                                               interval_intersect)
    end

    return nothing
end

"""
    forward_multiply_binary!

Updates storage tapes with forward evalution for node representing `n = x*y`.
"""
function forward_multiply_binary!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                                  numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                                  x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, subgrad_tol::Float64,
                                  sparsity::Vector{Int},
                                  is_post::Bool, is_intersect::Bool, is_first_eval::Bool, interval_intersect::Bool) where {N,T<:RelaxTag}
    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = one(MC{N,T})
        num1 = numberstorage[arg1_index]
    else
        num1 = 1.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index =  children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = one(MC{N,T})
    else
        set2 = setstorage[arg2_index]
        num2 = 1.0
    end

    output_is_number = arg1_is_number && arg2_is_number

    # a * b
    if output_is_number
         numberstorage[k] = num1 * num2

    # x * b
    elseif !arg1_is_number && arg2_is_number
        outset = set1*num2 #is_first_eval ? (set1 * num2) : mult_kernel(set1, num2, setstorage[k].Intv)

    # a * y
    elseif arg1_is_number && !arg2_is_number
        outset = num1*set2 #is_first_eval ? (num1 * set2) : mult_kernel(set2, num1, setstorage[k].Intv)

    # x * y
    else
        outset = set1*set2 #is_first_eval ? (set1 * set2) : mult_kernel(set1, set2, setstorage[k].Intv)

    end

    numvalued[k] = output_is_number
    if !output_is_number
         setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect,
                                                         interval_intersect)
    end

    return nothing
end


"""
    forward_multiply_narity!

Updates storage tapes with forward evalution of node representing `n = *(x, y, z,...)`.
"""
function forward_multiply_narity!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                                  numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                                  x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, subgrad_tol::Float64,
                                  sparsity::Vector{Int},
                                  is_post::Bool, is_intersect::Bool, interval_intersect::Bool) where {N,T<:RelaxTag}
    # get row indices
    idx = first(children_idx)

    # extract values for argument 1
    arg_index =  children_arr[idx]
    output_is_number =  numvalued[arg_index]
    if output_is_number
        tmp_set = 1.0#one(MC{N,T})
        tmp_num = numberstorage[arg_index]
    else
        tmp_num = 1.0
        tmp_set = setstorage[arg_index]
    end


    for idx = 2:length(children_idx)
        cidx =  children_idx[idx]
        arg_index_t =  children_arr[cidx]
        arg_is_number_t = numvalued[arg_index_t]
        if arg_is_number_t
            tmp_num = tmp_num*numberstorage[arg_index_t]
        else
            tmp_set = tmp_set*setstorage[arg_index_t]
        end
        output_is_number &= arg_is_number_t
    end

     numvalued[k] = output_is_number
    if output_is_number
         numberstorage[k] = tmp_num
    else
       tmp_set *= tmp_num
       setstorage[k] = overwrite_or_intersect(tmp_set, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect,
                                              interval_intersect)
    end

    return nothing
end

"""
    forward_minus!

Updates storage tapes with forward evalution for node representing `n = x-y`.
"""
function forward_minus!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                        numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                        x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, subgrad_tol::Float64,
                        sparsity::Vector{Int},
                        is_post::Bool, is_intersect::Bool, is_first_eval::Bool, interval_intersect::Bool) where {N,T<:RelaxTag}

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 =  setstorage[arg2_index]
        num2 = 0.0
    end

    output_is_number = arg1_is_number && arg2_is_number

    # a - b
    if output_is_number
         numberstorage[k] = num1 - num2

    # x - b
    elseif !arg1_is_number && arg2_is_number
        outset = set1 - num2 #is_first_eval ? (set1 - num2) : minus_kernel(set1, num2, setstorage[k].Intv)

    # a - y
    elseif arg1_is_number && !arg2_is_number
        outset = num1 - set2 #is_first_eval ? (num1 - set2) : minus_kernel(num1, set2, setstorage[k].Intv)

    # x - y
    else
        outset = set1 - set2 #is_first_eval ? (set1 - set2) : minus_kernel(set1, set2, setstorage[k].Intv)

    end

    numvalued[k] = output_is_number
    if !output_is_number
         setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect,
                                                         interval_intersect)
    end

    return nothing
end

"""
    forward_power!

Updates storage tapes with forward evalution for node representing `n = x^y`.
"""
function forward_power!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                        numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                        x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64},
                        subgrad_tol::Float64, sparsity::Vector{Int},
                        is_post::Bool, is_intersect::Bool, is_first_eval::Bool, interval_intersect::Bool,
                        ctx::GuardCtx) where {N,T<:RelaxTag}

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 =  numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 =  setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index =  children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = setstorage[arg2_index]
        num2 = 0.0
    end

    # is output a number (by closure of the reals)
    output_is_number = arg1_is_number && arg2_is_number
    numvalued[k] = output_is_number

    # x^1 = x
    if arg2_is_number && (num2 == 1.0)
        if arg1_is_number
             numberstorage[k] = num1
        else
             setstorage[k] = set1
        end
        return nothing

    # x^0 = 1
    elseif arg2_is_number && (num2 == 0.0)
        if arg1_is_number
             numberstorage[k] = 1.0
        else
             setstorage[k] = zero(MC{N,T})
        end
        return nothing

    else
        # a^b
        if arg1_is_number && arg2_is_number
             numberstorage[k] = num1^num2
        # x^b
        elseif !arg1_is_number && arg2_is_number
            outset = set1^num2
        # a^y
        elseif arg1_is_number  && !arg2_is_number
            outset = num1^set2 # overdub(ctx, pow, num1, set2)
        # x^y
        elseif !arg1_is_number && !arg2_is_number
            outset = set1^set2
        end
    end

    if !output_is_number
        setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect,
                                               interval_intersect)
    end

    return nothing
end

"""
    forward_divide!

Updates storage tapes with forward evalution for node representing `n = x/y`.
"""
function forward_divide!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                         numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                         x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, subgrad_tol::Float64,
                         sparsity::Vector{Int},
                         is_post::Bool, is_intersect::Bool, is_first_eval::Bool, interval_intersect::Bool) where {N,T<:RelaxTag}

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = setstorage[arg2_index]
        num2 = 0.0
    end

    # is output a number (by closure of the reals)?
    output_is_number = arg1_is_number && arg2_is_number
    numvalued[k] = output_is_number

    # a/b
    if output_is_number
        numberstorage[k] = num1/num2

    # x/b
    elseif !arg1_is_number && arg2_is_number
        outset = set1/num2

    # a/y
    elseif arg1_is_number && !arg2_is_number
        outset = num1/set2

    # x/y
    else
        outset = set1/set2
    end

    if !output_is_number
        setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, is_post,
                                               is_intersect, interval_intersect)
    end

    return nothing
end

"""
    forward_user_multivariate!

Updates storage tapes with forward evalution for node representing `n = user_f(x, y...)`.
"""
function forward_user_multivariate!(k::Int64, op::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                                    numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                                    x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, subgrad_tol::Float64,
                                    sparsity::Vector{Int},
                                    is_post::Bool, is_intersect::Bool, interval_intersect::Bool, ctx::GuardCtx,
                                    user_operators::JuMP._Derivatives.UserOperatorRegistry,
                                    num_mv_buffer::Vector{Float64}) where {N, T<:RelaxTag}

    n = length(children_idx)
    evaluator = user_operators.multivariate_operator_evaluator[op - JuMP._Derivatives.USER_OPERATOR_ID_START + 1]
    set_input = zeros(MC{N,T}, n)
    num_input = view(num_mv_buffer, 1:n)
    fill!(num_input, -Inf)

    buffer_count = 1
    output_is_number = true
    for c_idx in children_idx
        arg_index = children_arr[c_idx]
        arg_is_number = numvalued[arg_index]
        if arg_is_number
            num_input[buffer_count] = numberstorage[arg_index]
        else
            set_input[buffer_count] = setstorage[arg_index]
        end
        output_is_number &= arg_is_number
        buffer_count += 1
    end

    if output_is_number
        numberstorage[k] = MOI.eval_objective(evaluator, num_input)
    else
        for i = 1:(buffer_count - 1)
            if !isinf(num_input[i])
                 set_input[buffer_count] = MC{N,T}(num_input[buffer_count])
            end
        end
        outset = MOI.eval_objective(evaluator, set_input)
        setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect,
                                               interval_intersect)
    end
    numvalued[k] = output_is_number

    return nothing
end

"""
    forward_univariate_number!

Updates storage tapes with forward evalution for node representing `n = f(c)` where f is standard function
and `c` is a number.
"""
function forward_univariate_number!(k::Int64, op::Int64, arg_idx::Int, numvalued::Vector{Bool}, numberstorage::Vector{Float64})

    tmp_num = numberstorage[arg_idx]
    outnum = eval_univariate_set(op, tmp_num)

    numberstorage[k] = outnum
    numvalued[k] = true

    return nothing
end

"""
    forward_univariate_tiepnt_1!

Updates storage tapes with forward evalution for node representing `n = f(x)` where f is standard function
that requires a single tiepoint calculation per convex/concave relaxation (e.g. tan).
"""
function forward_univariate_tiepnt_1!(k::Int64, op::Int64, child_idx::Int64, setstorage::Vector{V},
                                      x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, subgrad_tol::Float64,
                                      sparsity::Vector{Int},
                                      tpdict::Dict{Int64, Tuple{Int64,Int64,Int64,Int64}},
                                      tp1storage::Vector{Float64}, tp2storage::Vector{Float64},
                                      is_post::Bool, is_intersect::Bool, is_first_eval::Bool, interval_intersect::Bool, ctx::GuardCtx) where V

    tmp_set = setstorage[child_idx]

    tidx1, tidx2 = tpdict[k]
    tp1 =  tp1storage[tidx1]
    tp2 =  tp2storage[tidx1]
    new_tie_points = tp1 === Inf

    outset, tp1, tp2 = single_tp_set(op, tmp_set, setstorage[k], tp1, tp2, is_first_eval)

    if new_tie_points
        tp1storage[tidx1] = tp1
        tp2storage[tidx1] = tp2
    end

    setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect,
                                           interval_intersect)
    return nothing
end

"""
    forward_univariate_tiepnt_2!

Updates storage tapes with forward evalution for node representing `n = f(x)` where f is standard function
that requires a two tiepoint calculations per convex/concave relaxation (e.g. sin).
"""
function forward_univariate_tiepnt_2!(k::Int64, op::Int64, child_idx::Int64, setstorage::Vector{V},
                                      x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64},
                                      subgrad_tol::Float64, sparsity::Vector{Int},
                                      tpdict::Dict{Int64, Tuple{Int64,Int64,Int64,Int64}},
                                      tp1storage::Vector{Float64}, tp2storage::Vector{Float64},
                                      tp3storage::Vector{Float64}, tp4storage::Vector{Float64},
                                      is_post::Bool, is_intersect::Bool, is_first_eval::Bool,
                                      interval_intersect::Bool, ctx::GuardCtx) where V

    tmp_set = setstorage[child_idx]

    # retreive previously calculated tie-points
    # These are re-initialize to Inf for each box
    tidx1, tidx2 = tpdict[k]
    tp1 = tp1storage[tidx1]
    tp2 = tp2storage[tidx1]
    tp3 = tp3storage[tidx2]
    tp4 = tp4storage[tidx2]

    new_tie_points = tp1 === Inf

    # Perform an evaluation of the univariate function overdubbed with Cassette.jl
    outset, tp1, tp2, tp3, tp4 = double_tp_set(op, tmp_set, setstorage[k], tp1, tp2, tp3, tp4, is_first_eval)

    # Store new tiepoints if new evaluation
    if new_tie_points
        tp1storage[tidx1] = tp1
        tp2storage[tidx1] = tp2
        tp3storage[tidx2] = tp3
        tp4storage[tidx2] = tp4
    end

    setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect,
                                           interval_intersect)
    return nothing
end

"""
    forward_univariate_user!

Updates storage tapes with forward evalution for node representing `n = user_f(x)`.
"""
function forward_univariate_user!(k::Int64, op::Int64, child_idx::Int64, arg_is_number::Bool,
                                  setstorage::Vector{V}, x::Vector{Float64}, lbd::Vector{Float64},
                                  ubd::Vector{Float64}, subgrad_tol::Float64, sparsity::Vector{Int},
                                  is_post::Bool, is_intersect::Bool, is_first_eval::Bool,
                                  interval_intersect::Bool, ctx::GuardCtx, user_operators) where V

    userop = op - JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START + 1
    f = user_operators.univariate_operator_f[userop]

    if arg_is_number
        tmp_num = setstorage[child_idx]
        outnum = f(tmp_num)
        numberstorage[k] = outnum

    else
        tmp_set = setstorage[child_idx]
        outset = f(tmp_set)
        setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity,
                                               is_post, is_intersect, interval_intersect)
    end

    return nothing
end

"""
    forward_univariate_other!

Updates storage tapes with forward evalution for node representing `n = f(x)` where `f` is a standard function
that does not require a tiepoint evaluation (e.g. exp).
"""
function forward_univariate_other!(k::Int64, op::Int64, child_idx::Int64, setstorage::Vector{V},
                                   x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64},
                                   subgrad_tol::Float64, sparsity::Vector{Int}, is_post::Bool,
                                   is_intersect::Bool, is_first_eval::Bool, interval_intersect::Bool,
                                   ctx::GuardCtx) where V

    tmp_set = setstorage[child_idx]
    outset = eval_univariate_set(op, tmp_set)
    setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x, lbd, ubd, subgrad_tol, sparsity,
                                           is_post, is_intersect, interval_intersect)

    return nothing
end

"""
    expand_set

"""
function expand_set(::Type{MC{N2,T}}, x::MC{N1,T}, fsparse::Vector{Int64},
                    subsparse::Vector{Int64}, cv_buffer::Vector{Float64},
                    cc_buffer::Vector{Float64}) where {N1 ,N2, T<:RelaxTag}

    # TODO: Perform this in a manner that does not allocate via an ntuple
    # operator or via generated code
    cvg = x.cv_grad
    ccg = x.cc_grad
    xcount = 1
    xcurrent = subsparse[1]
    for i = 1:N2
        if fsparse[i] === xcurrent
            cv_buffer[i] = cvg[xcount]
            cc_buffer[i] = ccg[xcount]
            xcount += 1
            if xcount <= N1
                xcurrent = subsparse[xcount]
            else
                break
            end
        else
            cv_buffer[i] = 0.0
        end
    end
    cv_grad = SVector{N2,Float64}(cv_buffer)
    cc_grad = SVector{N2,Float64}(cc_buffer)
    return MC{N2,T}(x.cv, x.cc, x.Intv, cv_grad, cc_grad, x.cnst)
end

"""
    forward_get_subexpression!

Unpacks the `MC{N1,T}` in the subexpression to a `MC{N2,T}` where `N1` is the sparsity of the
subexpression and `N2` is the sparsity of the function tape. Note that the sparsity of the
subexpression is less than the sparsity of the function itself. This limits the storage
required by tapes but prevents reverse mode subgradient propagation through expressions
[This may be subject to change in the future once the reverse mode propagation becomes
more robust].
"""
function forward_get_subexpression!(k::Int64, active_subexpr::NonlinearExpression{MC{N,T}},
                                    setstorage::Vector{MC{Q,T}}, numberstorage::Vector{Float64},
                                    numvalued::Vector{Bool}, fsparsity::Vector{Int64},
                                    cv_buffer::Vector{Float64}, cc_buffer::Vector{Float64}) where {N,Q,T<:RelaxTag}

    is_number = active_subexpr.isnumber[1]
    if is_number
        numberstorage[k] = active_subexpr.setstorage[1]
    else
        set_arg = active_subexpr.setstorage[1]
        setstorage[k] = expand_set(MC{Q,T}, set_arg, fsparsity, active_subexpr.grad_sparsity, cv_buffer, cc_buffer)
    end
    numvalued[k] = is_number

    return nothing
end

const id_to_operator = Dict(value => key for (key, value) in JuMP.univariate_operator_to_id)

"""
    forward_pass_kernel!

Performs a forward pass using the tape information passed as arguments. Each variety of node calls an associated
forward_xxx function where xxx is a descriptor.
"""
function forward_pass_kernel!(nd::Vector{JuMP.NodeData}, adj::SparseMatrixCSC{Bool,Int64}, x::Vector{Float64},
                              lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                              setstorage::Vector{MC{N,T}}, numberstorage::Vector{Float64}, numvalued::Vector{Bool},
                              tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}, tp1storage::Vector{Float64},
                              tp2storage::Vector{Float64}, tp3storage::Vector{Float64}, tp4storage::Vector{Float64},
                              user_operators::JuMP._Derivatives.UserOperatorRegistry, subexpressions::Vector{NonlinearExpression},
                              func_sparsity::Vector{Int64}, reverse_sparsity::Vector{Int64},
                              num_mv_buffer::Vector{Float64}, ctx::GuardCtx,
                              is_post::Bool, is_intersect::Bool, is_first_eval::Bool, interval_intersect::Bool,
                              cv_buffer::Vector{Float64}, cc_buffer::Vector{Float64},
                              treat_x_as_number::Vector{Bool}, subgrad_tol::Float64) where {N, T<:RelaxTag}

    children_arr = rowvals(adj)

    FORWARD_DEBUG && println(" ")
    for k = length(nd):-1:1

        oldset = setstorage[k]
        nod = nd[k]
        op = nod.index

        if nod.nodetype == JuMP._Derivatives.VALUE
            numvalued[k] = true
            FORWARD_DEBUG && println("value[$op]    at k = $k -> $(numberstorage[k])")

        elseif nod.nodetype == JuMP._Derivatives.PARAMETER
            numvalued[k] = true
            FORWARD_DEBUG && println("parameter[$op] at k = $k -> $(numberstorage[k])")

        elseif nod.nodetype == JuMP._Derivatives.VARIABLE
            isa_number =  treat_x_as_number[op]
            numvalued[k] = isa_number
            xval = x[op]
            if isa_number
                 numberstorage[k] = xval
            else
                seed_index = reverse_sparsity[op]
                seed_grad = seed_gradient(seed_index, Val{N}())
                lower_bnd_val = lbd[op]
                upper_bnd_val = ubd[op]
                xcv_eps = xval - (xval - lower_bnd_val)*1E-8
                xcc_eps = xval + (upper_bnd_val - xval)*1E-8
                xMC = MC{N,T}(xcv_eps, xcc_eps, Interval{Float64}(lower_bnd_val, upper_bnd_val),
                              seed_grad, seed_grad, false)
                setstorage[k] = is_first_eval ? xMC : (xMC ∩ setstorage[k].Intv)
            end
            FORWARD_DEBUG && println("variable[$op] at k = $k -> $(setstorage[k])")
        elseif nod.nodetype == JuMP._Derivatives.SUBEXPRESSION
            active_subexpr = subexpressions[op]
            forward_get_subexpression!(k, active_subexpr, setstorage, numberstorage, numvalued, func_sparsity,
                                       cv_buffer, cc_buffer)

        elseif nod.nodetype == JuMP._Derivatives.CALL

            children_idx = nzrange(adj, k)
            n_children = length(children_idx)

            # :+ with arity two or greater
            if op === 1
                n = length(children_idx)
                if n === 2
                    forward_plus_binary!(k, children_arr, children_idx, numvalued, numberstorage,
                                         setstorage, x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect, is_first_eval,
                                         interval_intersect)
                else
                    forward_plus_narity!(k, children_arr, children_idx, numvalued, numberstorage,
                                         setstorage, x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect, interval_intersect)
                end
                FORWARD_DEBUG && println("plus[$n]     at k = $k -> $(setstorage[k])")

            # :- with arity two
            elseif op === 2
                forward_minus!(k, children_arr, children_idx, numvalued, numberstorage,
                               setstorage, x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect, is_first_eval,
                               interval_intersect)
                FORWARD_DEBUG && println("minus        at k = $k -> $(setstorage[k])")

            # :* with arity two or greater
            elseif op === 3
                n = length(children_idx)
                if n === 2
                    forward_multiply_binary!(k, children_arr, children_idx, numvalued,
                                             numberstorage, setstorage, x, lbd, ubd, subgrad_tol, sparsity, is_post,
                                             is_intersect, is_first_eval, interval_intersect)
                else
                    forward_multiply_narity!(k, children_arr, children_idx, numvalued,
                                             numberstorage, setstorage, x, lbd, ubd, subgrad_tol,
                                             sparsity, is_post, is_intersect, interval_intersect)
                end
                FORWARD_DEBUG && println("mult[$n]     at k = $k -> $(setstorage[k])")

            # :^
            elseif op === 4
                forward_power!(k, children_arr, children_idx, numvalued, numberstorage,
                               setstorage, x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect, is_first_eval,
                               interval_intersect, ctx)

                FORWARD_DEBUG && println("power       at k = $k -> $(setstorage[k])")

            # :/
            elseif op === 5
                forward_divide!(k, children_arr, children_idx, numvalued, numberstorage,
                                setstorage, x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect, is_first_eval,
                                interval_intersect, ctx)

                FORWARD_DEBUG && println("divide      at k = $k -> $(setstorage[k])")

            # user multivariate function
            elseif op >= JuMP._Derivatives.USER_OPERATOR_ID_START
                forward_user_multivariate!(k, op, children_arr, children_idx, numvalued, numberstorage,
                                           setstorage, x, lbd, ubd, subgrad_tol, sparsity, is_post, is_intersect,
                                           interval_intersect, ctx, user_operators, num_mv_buffer)
                FORWARD_DEBUG && println("user_mult   at k = $k -> $(setstorage[k])")

            else
               error("Unsupported operation $(operators[op])")
            end

        elseif nod.nodetype == JuMP._Derivatives.CALLUNIVAR

            # checks to see if operator is a number
            child_idx = first(nzrange(adj, k))
            arg_idx = children_arr[adj.colptr[k]]
            arg_is_number =  numvalued[arg_idx]
            numvalued[k] = arg_is_number

            # performs univariate operators on number valued inputs
            if op >= JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START
                forward_univariate_user!(k, op, arg_idx, arg_is_number, setstorage, x, lbd, ubd, subgrad_tol,
                                         sparsity, is_post, is_intersect, is_first_eval, interval_intersect, ctx,
                                         user_operators)

            elseif arg_is_number
                forward_univariate_number!(k, op, arg_idx, numvalued, numberstorage)

            # performs set valued operators that require a single tiepoint calculation
            elseif single_tp(op)
                forward_univariate_tiepnt_1!(k, op, arg_idx, setstorage, x, lbd, ubd,  subgrad_tol, sparsity, tpdict,
                                             tp1storage, tp2storage, is_post, is_intersect,
                                             is_first_eval, interval_intersect, ctx)

            # performs set valued operators that require two tiepoint calculations
            elseif double_tp(op)
                forward_univariate_tiepnt_2!(k, op, arg_idx, setstorage, x, lbd, ubd,  subgrad_tol, sparsity, tpdict,
                                             tp1storage, tp2storage, tp3storage, tp4storage,
                                             is_post, is_intersect, is_first_eval,
                                             interval_intersect, ctx)

            # performs set valued operator on other functions in base library
            else
                forward_univariate_other!(k, op, arg_idx, setstorage, x, lbd, ubd, subgrad_tol, sparsity, is_post,
                                          is_intersect, is_first_eval, interval_intersect, ctx)

            end
            FORWARD_DEBUG && println("fop[$op]   at k = $k -> $(setstorage[k])")
        else
            error("Unrecognized node type $(nod.nodetype).")

        end

    end

    return nothing
end


###
### Define forward evaluation pass
###
function forward_pass!(evaluator::Evaluator, d::NonlinearExpression{V}) where V
    # check that prior subexpressions have been evaluated
    # i.e. box_id is same and reference point is the same
    for i = 1:d.dependent_subexpression_count
        if !prior_eval(evaluator, i)
            subexpr = evaluator.subexpressions[i]
            forward_pass!(evaluator, subexpr)
        end
    end

    forward_pass_kernel!(d.nd, d.adj, evaluator.x, evaluator.lower_variable_bounds,
                         evaluator.upper_variable_bounds, d.grad_sparsity,
                         d.setstorage,
                         d.numberstorage, d.isnumber, d.tpdict,
                         d.tp1storage, d.tp2storage, d.tp3storage, d.tp4storage,
                         evaluator.user_operators, evaluator.subexpressions,
                         d.grad_sparsity, d.reverse_sparsity,
                         evaluator.num_mv_buffer, evaluator.ctx,
                         evaluator.is_post, evaluator.is_intersect,
                         evaluator.is_first_eval, evaluator.interval_intersect,
                         d.cv_grad_buffer, d.cc_grad_buffer,
                         evaluator.treat_x_as_number, evaluator.subgrad_tol)
    return nothing
end

function forward_pass!(evaluator::Evaluator, d::BufferedNonlinearFunction{V}) where V

    forward_pass!(evaluator, d.expr)
    d.has_value = true
    d.last_past_reverse = false

    return nothing
end
