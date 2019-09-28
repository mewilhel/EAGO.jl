# Defines functions required for linear algebra packages
one(x::MC{N}) where N = MC{N}(1.0, 1.0, one(Interval{Float64}), zeros(SVector{N,Float64}), zeros(SVector{N,Float64}), true)
zero(x::MC{N}) where N = MC{N}(0.0, 0.0, zero(Interval{Float64}), zeros(SVector{N,Float64}), zeros(SVector{N,Float64}), true)
real(x::MC) = x
@inline dist(x1::MC, x2::MC) = max(abs(x1.cc - x2.cc), abs(x1.cv - x2.cv))
@inline eps(x::MC) = max(eps(x.cc), eps(x.cv))
@inline mid(x::MC) = mid(x.Intv)

function plus_kernel(x::MC{N}, y::MC{N}, z::Interval{Float64}) where N
	MC{N}(x.cv + y.cv, x.cc + y.cc, z, x.cv_grad + y.cv_grad, x.cc_grad + y.cc_grad, (x.cnst && y.cnst))
end
+(x::MC, y::MC) = plus_kernel(x, y, x.Intv + y.Intv)
plus_kernel(x::MC, y::Interval{Float64}) = x
+(x::MC) = x

minus_kernel(x::MC{N}, z::Interval{Float64}) where N = MC{N}(-x.cc, -x.cv, z, -x.cc_grad, -x.cv_grad, x.cnst)
-(x::MC) = minus_kernel(x, -x.Intv)

function minus_kernel(x::MC{N}, y::MC{N}, z::Interval{Float64}) where N
	MC{N}(x.cv - y.cc, x.cc - y.cv, z, x.cv_grad - y.cc_grad, x.cc_grad - y.cv_grad, (x.cnst && y.cnst))
end
-(x::MC{N}, y::MC{N}) where N = minus_kernel(x, y, x.Intv - y.Intv)

################## CONVERT THROUGH BINARY DEFINITIONS #########################

# Addition
function plus_kernel(x::MC{N}, y::Float64, z::Interval{Float64}) where N
	MC{N}(x.cv + y, x.cc + y, z, x.cv_grad, x.cc_grad, x.cnst)
end
function plus_kernel(y::Float64, x::MC{N}, z::Interval{Float64}) where N
	MC{N}(x.cv + y, x.cc + y, z, x.cv_grad, x.cc_grad, x.cnst)
end
+(x::MC, y::Float64) = plus_kernel(x, y, (x.Intv + y))
+(y::Float64, x::MC) = plus_kernel(x, y, (x.Intv + y))

plus_kernel(x::MC, y::C) where {C <: AbstractFloat} = plus_kernel(x, convert(Float64, y))
plus_kernel(x::C, y::MC) where {C <: AbstractFloat} = plus_kernel(convert(Float64, x), y)
+(x::MC, y::C) where {C <: AbstractFloat} = x + convert(Float64, y)
+(y::C, x::MC) where {C <: AbstractFloat} = x + convert(Float64, y)

plus_kernel(x::MC, y::C) where {C <: Integer} = plus_kernel(x, convert(Float64, y))
plus_kernel(x::C, y::MC) where {C <: Integer} = plus_kernel(convert(Float64, x), y)
+(x::MC, y::C) where {C <: Integer} = x + convert(Float64, y)
+(y::C, x::MC) where {C <: Integer} = x + convert(Float64, y)

# Subtraction
function minus_kernel(x::MC{N}, c::Float64, z::Interval{Float64}) where N
	MC{N}(x.cv - c, x.cc - c, z, x.cv_grad, x.cc_grad, x.cnst)
end
function minus_kernel(c::Float64, x::MC{N}, z::Interval{Float64}) where N
	MC{N}(c - x.cc, c - x.cv, z, -x.cc_grad, -x.cv_grad, x.cnst)
end
-(x::MC, c::Float64) = minus_kernel(x, c, x.Intv - c)
-(c::Float64, x::MC) = minus_kernel(c, x, c - x.Intv)

minus_kernel(x::MC, y::C, z::Interval{Float64}) where {C <: AbstractFloat} = minus_kernel(x, convert(Float64,y), z)
minus_kernel(y::C, x::MC, z::Interval{Float64}) where {C <: AbstractFloat} = minus_kernel(convert(Float64,y), x, z)
-(x::MC, c::C) where {C <: AbstractFloat} = minus_kernel(x, convert(Float64,c), x.Intv)
-(c::C, x::MC) where {C <: AbstractFloat} = minus_kernel(convert(Float64,c), x, x.Intv)

minus_kernel(x::MC, y::C, z::Interval{Float64}) where {C <: Integer} = minus_kernel(x, convert(Float64,y), z)
minus_kernel(y::C, x::MC, z::Interval{Float64}) where {C <: Integer} = minus_kernel(convert(Float64,y), x, z)
-(x::MC, c::C) where {C <: Integer} = minus_kernel(x, convert(Float64,c), x.Intv)
-(c::C, x::MC) where {C <: Integer} = minus_kernel(convert(Float64,c), x, x.Intv)

# Multiplication
function mult_kernel(x::MC{N}, c::Float64, z::Interval{Float64}) where N
	if (c >= 0.0)
		z = MC{N}(c*x.cv, c*x.cc, z, c*x.cv_grad, c*x.cc_grad, x.cnst)
	else
		z = MC{N}(c*x.cc, c*x.cv, z, c*x.cc_grad, c*x.cv_grad, x.cnst)
	end
	return z
end
mult_kernel(c::Float64, x::MC, z::Interval{Float64}) = mult_kernel(x, c, z)
*(x::MC, c::Float64) = mult_kernel(x, c, c*x.Intv)
*(c::Float64, x::MC) = mult_kernel(x, c, c*x.Intv)

mult_kernel(x::MC, c::C, z::Interval{Float64}) where {C<:AbstractFloat} = mult_kernel(x, convert(Float64, c), z)
mult_kernel(c::C, x::MC, z::Interval{Float64}) where {C<:AbstractFloat} =  mult_kernel(x, convert(Float64, c), z)
*(c::C, x::MC) where {C<:AbstractFloat}  = x*Float64(c)
*(x::MC, c::C) where {C<:AbstractFloat}  = x*Float64(c)

mult_kernel(x::MC, c::C, z::Interval{Float64}) where {C<:Integer} = mult_kernel(x, convert(Float64, c), z)
mult_kernel(c::C, x::MC, z::Interval{Float64}) where {C<:Integer} =  mult_kernel(x, convert(Float64, c), z)
*(c::C, x::MC) where {C<:Integer}  = x*Float64(c)
*(x::MC, c::C) where {C<:Integer}  = x*Float64(c)

# Division
div_kernel(x::MC, y::C, z::Interval{Float64}) where {C<:AbstractFloat}  = mult_kernel(x, inv(y), z)
div_kernel(x::C, y::MC, z::Interval{Float64}) where {C<:AbstractFloat}  = mult_kernel(x, inv(y), z)
div_kernel(x::MC, y::C, z::Interval{Float64}) where {C<:Integer}  = mult_kernel(x, inv(y), z)
div_kernel(x::C, y::MC, z::Interval{Float64}) where {C<:Integer}  = mult_kernel(x, inv(y), z)
/(x::MC, y::C) where {C<:AbstractFloat} = x*inv(convert(Float64,y))
/(x::C, y::MC) where {C<:AbstractFloat} = convert(Float64,x)*inv(y)
/(x::MC, y::C) where {C<:Integer} = x*inv(convert(Float64,y))
/(x::C, y::MC) where {C<:Integer} = convert(Float64,x)*inv(y)

# Maximization
max_kernel(c::Float64, x::MC, z::Interval{Float64}) = max_kernel(x, c, z)
max_kernel(x::MC, c::C, z::Interval{Float64}) where {C<:AbstractFloat} = max_kernel(x, convert(Float64, c), z)
max_kernel(c::C, x::MC, z::Interval{Float64}) where {C<:AbstractFloat} = max_kernel(x, convert(Float64, c), z)
max_kernel(x::MC, c::C, z::Interval{Float64}) where {C<:Integer} = max_kernel(x, convert(Float64, c), z)
max_kernel(c::C, x::MC, z::Interval{Float64}) where {C<:Integer} = max_kernel(x, convert(Float64, c), z)

max(c::Float64, x::MC) = max_kernel(x, c, max(x.Intv, c))
max(x::MC, c::C) where {C<:AbstractFloat} = max_kernel(x, convert(Float64, c), max(x.Intv, c))
max(c::C, x::MC) where {C<:AbstractFloat} = max_kernel(x, convert(Float64, c), max(x.Intv, c))
max(x::MC, c::C) where {C<:Integer} = max_kernel(x, convert(Float64, c), max(x.Intv, c))
max(c::C, x::MC) where {C<:Integer} = max_kernel(x, convert(Float64, c), max(x.Intv, c))

# Minimization
min_kernel(x::MC, c::C, z::Interval{Float64}) where {C<:AbstractFloat} = min_kernel(x, convert(Float64, c), z)
min_kernel(c::C, x::MC, z::Interval{Float64}) where {C<:AbstractFloat} = min_kernel(x, convert(Float64, c), z)
min_kernel(x::MC, c::C, z::Interval{Float64}) where {C<:Integer} = min_kernel(x, convert(Float64, c), z)
min_kernel(c::C, x::MC, z::Interval{Float64}) where {C<:Integer} = min_kernel(x, convert(Float64, c), z)

min(c::Float64, x::MC) = min_kernel(x, c, min(x.Intv, c))
min(x::MC, c::C) where {C<:AbstractFloat} = min_kernel(x, convert(Float64, c), min(x.Intv, c))
min(c::C, x::MC) where {C<:AbstractFloat} = min_kernel(x, convert(Float64, c), min(x.Intv, c))
min(x::MC, c::C) where {C<:Integer} = min_kernel(x, convert(Float64, c), min(x.Intv, c))
min(c::C, x::MC) where {C<:Integer} = min_kernel(x, convert(Float64, c), min(x.Intv, c))

# Promote and Convert
promote_rule(::Type{MC{N}}, ::Type{S}) where {S<:Integer,N} = MC{N}
promote_rule(::Type{MC{N}}, ::Type{S}) where {S<:AbstractFloat,N} = MC{N}
promote_rule(::Type{MC{N}}, ::Type{S}) where {S<:Interval,N} = MC{N}
promote_rule(::Type{MC{N}}, ::Type{S}) where {S<:Real,N} = MC{N}

convert(::Type{MC{N}}, x::S) where {S<:Integer, N} = MC{N}(Interval{Float64}(x))
convert(::Type{MC{N}}, x::S) where {S<:AbstractFloat, N} = MC{N}(Interval{Float64}(x))
convert(::Type{MC{N}}, x::S) where {S<:Interval, N} = MC{N}(Interval{Float64}(x.lo, x.hi))
