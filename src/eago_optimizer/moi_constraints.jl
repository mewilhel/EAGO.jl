# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/constraints.jl
# Defines constraints supported Optimizer and how to store them.
#############################################################################


##### Supports function and add_constraint for scalar affine functions
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_linear(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m._input_problem.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_linear SAF LT _linear_leq_constraints _linear_leq_count
@define_addconstraint_linear SAF GT _linear_geq_constraints _linear_geq_count
@define_addconstraint_linear SAF ET _linear_eq_constraints _linear_eq_count

##### Supports function and add_constraint for scalar quadratic functions
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_quadratic(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m._input_problem.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_quadratic SQF LT _quadratic_leq_constraints _quadratic_leq_count
@define_addconstraint_quadratic SQF GT _quadratic_geq_constraints _quadratic_geq_count
@define_addconstraint_quadratic SQF ET _quadratic_eq_constraints _quadratic_eq_count

##### Supports function and add_constraint for conic functions
const CONE_SETS = Union{SECOND_ORDER_CONE, POWER_CONE, EXP_CONE}
MOI.supports_constraint(::Optimizer, ::Type{VECOFVAR}, ::Type{S}) where {S <: CONE_SETS} = true

macro @define_addconstraint_conic(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            if length(func.variables) !== set.dimension
                error("Dimension of $(s) does not match number of terms in $(f)")
            end
            check_inbounds!(m, func)
            push!(m._input_problem.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_quadratic VECOFVAR SECOND_ORDER_CONE _conic_second_order _conic_second_order_count
@define_addconstraint_quadratic VECOFVAR POWER_CONE        _conic_power        _conic_power_count
@define_addconstraint_quadratic VECOFVAR EXP_CONE          _conic_exp          _conic_exp_count
