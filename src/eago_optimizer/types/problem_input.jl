"""
$(TYPEDEF)

A structure used to hold objectives and constraints added to EAGO model.
The constraints generally aren't used for relaxations. Essentially, this structure
exists to cache a previous model so that alternative representations of subproblems
may then be constructed and assessed (i.e. performing a local solve in full-space
while adding auxillary variable to form a relaxation).
"""
Base.@kwdef mutable struct InputProblem <: MOI.ModelLike

    # TODO: Replace Dict with MOIU.VectorOfConstraint when version 0.9.21 of MOI is tagged
    _variable_leq::Dict{CI{SV,LT},Tuple{SV,LT}} = Dict{CI{SV,LT},Tuple{SV,LT}}()
    _variable_geq::Dict{CI{SV,GT},Tuple{SV,GT}} = Dict{CI{SV,GT},Tuple{SV,GT}}()
    _variable_eq::Dict{CI{SV,ET},Tuple{SV,ET}} = Dict{CI{SV,ET},Tuple{SV,ET}}()
    _variable_zo::Dict{CI{SV,ZO},Tuple{SV,ZO}} = Dict{CI{SV,ZO},Tuple{SV,ZO}}()
    _variable_int::Dict{CI{SV,MOI.Integer},Tuple{SV,MOI.Integer}} = Dict{CI{SV,MOI.Integer},Tuple{SV,MOI.Integer}}()

    # linear constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _linear_leq::Dict{CI{SAF,LT},Tuple{SAF,LT}} = Dict{CI{SAF,LT}, Tuple{SAF,LT}}()
    _linear_geq::Dict{CI{SAF,GT},Tuple{SAF,GT}} = Dict{CI{SAF,GT}, Tuple{SAF,GT}}()
    _linear_eq::Dict{CI{SAF,ET},Tuple{SAF,ET}} = Dict{CI{SAF,ET}, Tuple{SAF,ET}}()

    # quadratic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _quadratic_leq::Dict{CI{SQF,LT},Tuple{SQF,LT}} = Dict{CI{SQF,LT},Tuple{SQF,LT}}()
    _quadratic_geq::Dict{CI{SQF,GT},Tuple{SQF,GT}} = Dict{CI{SQF,GT},Tuple{SQF,GT}}()
    _quadratic_eq::Dict{CI{SQF,ET},Tuple{SQF,ET}} = Dict{CI{SQF,ET},Tuple{SQF,ET}}()

    # conic constraint storage and count (set by MOI.add_constraint in moi_constraints.jl)
    _conic_socp::Dict{CI{VECVAR,SOC_CONE},Tuple{VECVAR,SOC_CONE}} = Dict{CI{VECVAR, SOC_CONE}, Tuple{VECVAR, SOC_CONE}}()
    _conic_sdp::Dict{CI{VECVAR,PSD_CONE},Tuple{VECVAR,PSD_CONE}} = Dict{CI{VECVAR, PSD_CONE},Tuple{VECVAR, PSD_CONE}}()
    # nonlinear constraint storage
    _nonlinear_count::Int             = 0

    _variable_num::Int64                 = 0
    _constraint_index_num::Int           = 0
    _constraint_row_num::Int             = 0
    _constraint_offset::Vector{Int}      = Int[]

    # objective information (set by MOI.set(m, ::ObjectiveFunction...) in optimizer.jl)
    _objective_sv::SV = SV(VI(-1))
    _objective_saf::SAF = SAF(SAT[], 0.0)
    _objective_sqf::SQF = SQF(SAT[], SQT[], 0.0)
    _objective_type::ObjectiveType = UNSET

    # nlp constraints (set by MOI.set(m, ::NLPBlockData...) in optimizer.jl)
    _nlp_data::Union{MOI.NLPBlockData,Nothing} = nothing

    # objective sense information (set by MOI.set(m, ::ObjectiveSense...) in optimizer.jl)
    _optimization_sense::MOI.OptimizationSense = MOI.FEASIBILITY_SENSE
end

MOI.supports_constraint(::InputProblem,
                        ::Type{<:Union{SV, SAF, SQF}},
                        ::Type{<:Union{ET, GT, LT}},
                        ) = true

MOI.supports_constraint(::InputProblem,
                        ::Type{SV},
                        ::Type{<:Union{ZO, MOI.Integer}}) = true

MOI.supports_constraint(::InputProblem,
                        ::Type{<:Union{VECVAR}},
                        ::Type{<:Union{SOC_CONE, PSD_CONE}},
                        ) = true

MOI.supports(::InputProblem,
             ::Union{MOI.ObjectiveSense,
                     MOI.NLPBlock,
                     MOI.ObjectiveFunctionType,
                     MOI.ObjectiveFunction{SV},
                     MOI.ObjectiveFunction{SAF},
                     MOI.ObjectiveFunction{SQF}},
                     ) = true


@inline _variable_num(d::InputProblem) = d._variable_num
@inline function _integer_variable_num(d::InputProblem)
    intvars = first.(values(d._variable_zo))
    append!(intvars, first.(values(d._variable_int)))
    unique!(intvars)
    return length(intvars)
end
@inline function _second_order_cone_num(d::InputProblem)
    return length(d._conic_socp)
end
@inline function _psd_cone_num(d::InputProblem)
    return length(d._conic_sdp)
end
@inline function _quadratic_num(d::InputProblem)
    num = length(d._quadratic_leq)
    num += length(d._quadratic_geq)
    num += length(d._quadratic_eq)
    return num
end
@inline function _nl_expr_num(d::InputProblem)
    num = d._objective_type === NONLINEAR ? 1 : 0
    num += d._nonlinear_count
    return num
end
@inline _optimization_sense(d::InputProblem) = d._optimization_sense
@inline _objective_type(d::InputProblem)     = d._objective_type

function _check_inbounds!(d::InputProblem, val::T) where T <: Number
    if !(1 <= val <= d._variable_num)
        error("Invalid variable index $vi. ($(m._input_problem._variable_num) variables in the model.)")
    end
    return nothing
end
_check_inbounds!(d::InputProblem, vi::VI) = _check_inbounds!(d, vi.value)
_check_inbounds!(d::InputProblem, var::SV) = _check_inbounds!(d, var.variable)
_check_inbounds!(d::InputProblem, aff::SAF) = foreach(x -> _check_inbounds!(d, x.variable_index), aff.terms)
function _check_inbounds!(d::InputProblem, quad::SQF)
    foreach(x -> _check_inbounds!(d, x.variable_index), quad.affine_terms)
    for term in quad.quadratic_terms
        _check_inbounds!(d, term.variable_index_1)
        _check_inbounds!(d, term.variable_index_2)
    end
    return nothing
end
_check_inbounds!(d::InputProblem, v::VECVAR) = foreach(x -> _check_inbounds!(d, x), v.variables)

function MOI.add_variable(d::InputProblem)
    d._variable_num += 1
    return VI(d._variable_num)
end

const VARIABLE_TO_DICT = ((ET, :_variable_eq),  (LT, :_variable_leq),
                          (GT, :_variable_geq), (ZO, :_variable_zo),
                          (MOI.Integer, :_variable_int))

for (S, storage) in VARIABLE_TO_DICT
    for (T, Tdict) in VARIABLE_TO_DICT
        if storage == T
            @eval function MOI.set(d::InputProblem, ::MOI.ConstraintFunction, ci::CI{SV,$S}, f::$S)
                d.$(storage)[ci] = (f, d.$(storage)[ci][2])
                return
            end
            @eval function MOI.set(d::InputProblem, ::MOI.ConstraintSet, ci::CI{SV,$S}, s::$S)
                d.$(storage)[ci] = (d.$(storage)[ci][1], s)
                return
            end
        else
            @eval function MOI.set(d::InputProblem, ::MOI.ConstraintFunction, ci::CI{SV,$S}, f::$T)
                s = d.$(storage)[ci][2]
                delete!(d.$(storage), ci)
                d.$(Tdict)[ci] = (f, s)
                return
            end
            @eval function MOI.set(d::InputProblem, ::MOI.ConstraintSet, ci::CI{SV,$S}, s::$T)
                f = d.$(storage)[ci][1]
                delete!(d.$(storage), ci)
                d.$(Tdict)[ci] = (f, s)
                return
            end
        end
    end
    @eval function MOI.get(d::InputProblem, ::MOI.ListOfConstraintIndices{SV,$S})
        return collect(keys(d.$storage))
    end
    @eval function MOI.get(d::InputProblem, ::MOI.ConstraintFunction, ci::CI{SV,$S})
        return d.$storage[ci][1]
    end
    @eval function MOI.get(d::InputProblem, ::MOI.ConstraintSet, ci::CI{SV,$S})
        return d.$storage[ci][2]
    end
    @eval function MOI.get(d::InputProblem, ::MOI.ListOfConstraintAttributesSet{SV,$S})
        return MOI.AbstractConstraintAttribute[]
    end
    @eval function MOI.add_constraint(d::InputProblem, v::SV, s::$S)
        vi = v.variable
        _check_inbounds!(d, vi)
        ci = CI{SV, $S}(vi.value)
        d.$storage[ci] = (v, s)
        return ci
    end
end

macro define_constraint(F, S, dict_name)
    esc(quote
            function MOI.add_constraint(d::InputProblem, f::$F, s::$S)
                _check_inbounds!(d, f)
                d._constraint_index_num += 1
                ci = CI{$F,$S}(d._constraint_index_num)
                d.$(dict_name)[ci] = (copy(f), s)
                d._constraint_row_num += MOI.dimension(s)
                push!(d._constraint_offset, d._constraint_row_num)
                return ci
            end
            function MOI.get(d::InputProblem, ::MOI.ListOfConstraintIndices{$F,$S})
                return collect(keys(d.$dict_name))
            end
            function ($dict_name)(d::InputProblem)
                return d.$(dict_name)
            end
            function MOI.get(d::InputProblem, ::MOI.ConstraintFunction, ci::CI{$F,$S})
                return d.$dict_name[ci][1]
            end
            function MOI.get(d::InputProblem, ::MOI.ConstraintSet, ci::CI{$F,$S})
                return d.$dict_name[ci][2]
            end
            function MOI.get(d::InputProblem, ::MOI.ListOfConstraintAttributesSet{$F,$S})
                return MOI.AbstractConstraintAttribute[]
            end
        end)
end

macro define_objective(F, field_name)
    esc(quote
            function MOI.set(d::InputProblem, ::MOI.ObjectiveFunction{$F}, f::$F)
                _check_inbounds!(d, f)
                d.$field_name = copy(f)
                d._objective_type = _moi_to_obj_type(f)
                return nothing
            end
            function MOI.get(d::InputProblem, ::MOI.ObjectiveFunction{$F})
                return d.$field_name
            end
    end)
end

@define_constraint SAF LT _linear_leq
@define_constraint SAF GT _linear_geq
@define_constraint SAF ET _linear_eq

@define_constraint SQF LT _quadratic_leq
@define_constraint SQF GT _quadratic_geq
@define_constraint SQF ET _quadratic_eq
@define_constraint VECVAR SOC_CONE _conic_socp
@define_constraint VECVAR PSD_CONE _conic_sdp

_moi_to_obj_type(d::SV) = SINGLE_VARIABLE
_moi_to_obj_type(d::SAF) = SCALAR_AFFINE
_moi_to_obj_type(d::SQF) = SCALAR_QUADRATIC
