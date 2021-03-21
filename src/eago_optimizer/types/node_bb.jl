# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/node_bb.jl
# Defines storage for a node in the B&B tree & utilities functions
#############################################################################

"""
$(TYPEDEF)

Stores information associated with each node in Branch & Bound tree.

$(TYPEDFIELDS)
"""
struct NodeBB{N,T<:AbstractFloat}
    "Lower bounds of variable box."
    lower_variable_bound::NTuple{N,T}
    "Upper bounds of variable box."
    upper_variable_bound::NTuple{N,T}
    "Lower bound of problem solution on nodeBB"
    lower_bound::T
    "Upper bound of problem solution on nodeBB"
    upper_bound::T
    "Depth of node in B&B tree."
    depth::Int
    "Unique id for each node."
    id::Int
end

# Constructors
NodeBB{N,T}() = NodeBB(ntuple(x->zero(T), N), ntuple(x->zero(T), N), -Inf, Inf, 0, 1)
NodeBB(x::NodeBB{N,T}) = NodeBB{N,T}(x.lower_variable_bound, x.upper_variable_bound,
                                     x.lower_bound, x.upper_bound, x.depth, x.id)

# Copy utilities
Base.copy(x::NodeBB) = NodeBB(x.lower_variable_bound,
                              x.upper_variable_bound,
                              x.lower_bound,
                              x.upper_bound,
                              x.depth,
                              x.id)

# Access functions for broadcasting data easily
lower_variable_bound(x::NodeBB) = x.lower_variable_bound
upper_variable_bound(x::NodeBB) = x.upper_variable_bound
lower_variable_bound(x::NodeBB, i::Int) = x.lower_variable_bound[i]
upper_variable_bound(x::NodeBB, i::Int) = x.upper_variable_bound[i]
lower_variable_bound(x::NodeBB, id::Int, nid::Int) = x.lower_variable_bound[id:nid]
upper_variable_bound(x::NodeBB, id::Int, nid::Int) = x.upper_variable_bound[id:nid]
lower_bound(x::NodeBB) = x.lower_bound
upper_bound(x::NodeBB) = x.upper_bound
depth(x::NodeBB) = x.depth

# Iterations Functions
Base.isless(x::NodeBB, y::NodeBB) = x.lower_bound < y.lower_bound
Base.length(x::NodeBB{N,T}) where {N,T} = N
function Base.isempty(x::NodeBB{N,T}) where {N,T}
    for i = 1:N
        lower = @inbounds x.lower_variable_bound[i]
        upper = @inbounds x.upper_variable_bound[i]
        (lower > upper) && (return true)
    end
    return false
end

"""
$(FUNCTIONNAME)

Checks that node `x` and `y` have equal domains withing a tolerance of `atol`.
"""
function same_box(x::NodeBB{N,T}, y::NodeBB{N,T}, r::Float64) where {N,T}
    for i = 1:N
        ~isapprox(x.lower_variable_bound[i], y.lower_variable_bound[i], atol=r) && (return false)
        ~isapprox(x.upper_variable_bound[i], y.upper_variable_bound[i], atol=r) && (return false)
    end
    return true
end

# Compute middle & diameter
diam(x::NodeBB) = x.upper_variable_bound - x.lower_variable_bound
mid(x::NodeBB) = 0.5*(x.upper_variable_bound + x.lower_variable_bound)
mid(x::NodeBB, i) = 0.5*(x.upper_variable_bound[i] + x.lower_variable_bound[i])
