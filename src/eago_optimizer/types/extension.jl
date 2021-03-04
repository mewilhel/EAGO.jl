"""
$(TYPEDEF)

An abstract type the subtypes of which are associated with functions method
overloaded for for new extensions. An instance of the `DefaultExt <:ExtensionType`
structure to the `Optimizer` in the `ext_type` field.
"""
abstract type ExtensionType end
struct DefaultExt <: ExtensionType end
MOIU.map_indices(::Function, x::ExtensionType) = x
MOIU.map_indices(::Function, x::DefaultExt) = x

@enum(ObjectiveType, UNSET, SINGLE_VARIABLE, SCALAR_AFFINE, SCALAR_QUADRATIC, NONLINEAR)
@enum(ProblemType, UNCLASSIFIED, LP, MILP, SOCP, MISOCP, DIFF_CVX, MINCVX)
