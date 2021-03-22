mutable struct LinearEliminationMap
    a::Vector{Float64}                        = Float64[]
    b::Vector{Float64}                        = Float64[]
    c::Vector{Float64}                        = Float64[]
    xi::Vector{MOI.VariableIndex}             = MOI.VariableIndex[]
    xj::Vector{MOI.VariableIndex}             = MOI.VariableIndex[]
    eliminated::OrderedDict{CI{SAF,ET}, Bool} = Tuple{SAF,ET}()
end

UnitUpperTriangular

Ax = b -> (UU)x = b*
