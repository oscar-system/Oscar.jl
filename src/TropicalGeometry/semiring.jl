################################################################################
#
#  Tropical semirings and tropical semiring elements (= tropical numbers)
#
################################################################################

# minOrMax distinguishes between the min-plus and max-plus semiring
struct TropicalSemiring{minOrMax<:Union{typeof(min),typeof(max)}} <: Field
end

struct TropicalSemiringElem{minOrMax<:Union{typeof(min),typeof(max)}} <: FieldElem
    parent::TropicalSemiring{minOrMax}
    isinf::Bool        # distinguishes between ±∞ and other tropical numbers
    data::QQFieldElem

    function TropicalSemiringElem(R::TropicalSemiring{minOrMax}, isinf::Bool) where {minOrMax<:Union{typeof(min),typeof(max)}}
        @assert isinf
        return new{minOrMax}(R, true)
    end

    function TropicalSemiringElem(R::TropicalSemiring{minOrMax}, x::RingElem) where {minOrMax<:Union{typeof(min),typeof(max)}}
        return new{minOrMax}(R, false, x)
    end
end

# Type gymnastics
elem_type(::Type{TropicalSemiring{minOrMax}}) where {minOrMax<:Union{typeof(min),typeof(max)}} = TropicalSemiringElem{minOrMax}

parent_type(::Type{TropicalSemiringElem{minOrMax}}) where {minOrMax<:Union{typeof(min),typeof(max)}} = TropicalSemiring{minOrMax}



################################################################################
#
#  Basic access for tropical numbers
#
################################################################################

parent(x::TropicalSemiringElem) = x.parent

data(x::TropicalSemiringElem) = x.data    # undefined if x==±∞

isinf(x::TropicalSemiringElem) = x.isinf  # also serves as test if x==±∞

@doc raw"""
    convention(T::TropicalSemiring)

Return `min` if `T` is the min tropical semiring,
return `max` if `T` is the max tropical semiring.
Works similarly for tropical numbers,
tropical vectors and matrices, and tropical polynomials.

# Examples
```jldoctest; filter = r"\(generic function with .* methods\)"
julia> T = tropical_semiring(min)
Min tropical semiring

julia> convention(T)
min (generic function with 27 methods)

julia> T = tropical_semiring(max)
Max tropical semiring

julia> convention(T)
max (generic function with 27 methods)
```
"""
convention(T::TropicalSemiring{typeof(min)}) = min
convention(T::TropicalSemiring{typeof(max)}) = max
convention(a::TropicalSemiringElem{typeof(min)}) = min
convention(a::TropicalSemiringElem{typeof(max)}) = max
convention(v::Vector{TropicalSemiringElem{typeof(min)}}) = min
convention(v::Vector{TropicalSemiringElem{typeof(max)}}) = max
convention(M::Vector{Vector{TropicalSemiringElem{typeof(min)}}}) = min
convention(M::Vector{Vector{TropicalSemiringElem{typeof(max)}}}) = max
convention(M::Generic.MatSpaceElem{TropicalSemiringElem{typeof(min)}}) = min
convention(M::Generic.MatSpaceElem{TropicalSemiringElem{typeof(max)}}) = max
convention(f::Generic.MPoly{TropicalSemiringElem{typeof(min)}}) = min
convention(f::Generic.MPoly{TropicalSemiringElem{typeof(max)}}) = max


################################################################################
#
#  Constructors for tropical semirings
#
################################################################################

@doc raw"""
    tropical_semiring(M::Union{typeof(min),typeof(max)}=min)

Return the min-plus (default) or max-plus semiring.

!!! warning
    - `+`, `*`, `/`, and `^` are used for tropical addition, tropical multipliciation, tropical division, and tropical exponentiation, respectively.
    - There is no additive inverse or subtraction in the tropical semiring. Negating a tropical number or subtracting two tropical numbers will raise an error.
    - Zeroes of tropical semirings are printed as `infty` or `-infty` instead of their proper unicode characters.  To enabled unicode in the current and future sessions, run `allow_unicode(true)`.

# Examples (basic arithmetic)
```jldoctest
julia> T = tropical_semiring() # = tropical_semiring(min)
Min tropical semiring

julia> T = tropical_semiring(max)
Max tropical semiring

julia> 0*T(3) + 1*T(1)^2 + zero(T) # = max(0+3,1+2*1,-∞)
(3)

julia> T(0) == 0    # checks whether the tropical number is 0
true

julia> iszero(T(0)) # checks whether the tropical number is neutral element of addition
false
```

# Examples (polynomials)
```jldoctest
julia> T = tropical_semiring()
Min tropical semiring

julia> Tx,(x1,x2) = polynomial_ring(T,2)
(Multivariate polynomial ring in 2 variables over min tropical semiring, AbstractAlgebra.Generic.MPoly{TropicalSemiringElem{typeof(min)}}[x1, x2])

julia> f = x1 + -1*x2 + 0
x1 + (-1)*x2 + (0)

julia> evaluate(f,T.([-1//2,1//2])) # warning: omitting T() gives an error
(-1//2)
```

# Examples (matrices)
```jldoctest
julia> T = tropical_semiring()
Min tropical semiring

julia> A = identity_matrix(T, 2) # = tropical identity matrix
[  (0)   infty]
[infty     (0)]

julia> 2*A
[  (2)   infty]
[infty     (2)]

julia> A*A
[  (0)   infty]
[infty     (0)]

julia> det(A)
(0)

julia> minors(A,1)
4-element Vector{TropicalSemiringElem{typeof(min)}}:
 (0)
 infty
 infty
 (0)
```
"""
tropical_semiring() = TropicalSemiring{typeof(min)}()
tropical_semiring(::typeof(max)) = TropicalSemiring{typeof(max)}()
tropical_semiring(::typeof(min)) = TropicalSemiring{typeof(min)}()



################################################################################
#
#  Constructors for tropical numbers
#
################################################################################

function (T::TropicalSemiring)(u::TropicalSemiringElem)
  @req parent(u)==T "incompatible conventions"
  return u
end

zero(T::TropicalSemiring) = TropicalSemiringElem(T, true)    # neutral element w.r.t. +

one(T::TropicalSemiring) = TropicalSemiringElem(T, zero(QQ)) # neutral element w.r.t. *

inf(T::TropicalSemiring) = zero(T)

(T::TropicalSemiring)() = zero(T)



################################################################################
#
#  Conversion between tropical numbers and rational numbers.
#  If preserve_ordering==true and minOrMax==typeof(max), flip signs
#  (for info on tropical semiring ordering see comparison below).
#
################################################################################

function (R::TropicalSemiring{typeof(min)})(x::QQFieldElem; preserve_ordering::Bool=false)
  return TropicalSemiringElem(R,x)
end

function (R::TropicalSemiring{typeof(max)})(x::QQFieldElem; preserve_ordering::Bool=false)
  return (preserve_ordering ? TropicalSemiringElem(R,-x) : TropicalSemiringElem(R,x))
end

function (R::TropicalSemiring)(x::Union{Integer, Rational}; preserve_ordering::Bool=false)
  return R(QQ(x),preserve_ordering=preserve_ordering)
end

function (R::TropicalSemiring)(x::RingElem; preserve_ordering::Bool=false)
  x = QQ(x)
  @req parent(x)==QQ "cannot convert object of type $(repr(typeof(x)))"
  return R(x, preserve_ordering=preserve_ordering)
end

function (::QQField)(x::TropicalSemiringElem{typeof(min)}; preserve_ordering::Bool=false)
    @req !iszero(x) "cannot convert $(repr(x))"
    return data(x)
end

function (::QQField)(x::TropicalSemiringElem{typeof(max)}; preserve_ordering::Bool=false)
    @req !iszero(x) "cannot convert $(repr(x))"
    return (preserve_ordering ? -data(x) : data(x))
end

function (::ZZRing)(x::TropicalSemiringElem; preserve_ordering::Bool=false)
    return ZZ(QQ(x,preserve_ordering=preserve_ordering))
end

function (::Type{Int})(x::TropicalSemiringElem; preserve_ordering::Bool=false)
    return Int(QQ(x,preserve_ordering=preserve_ordering))
end

function (::Type{Rational})(x::TropicalSemiringElem; preserve_ordering::Bool=false)
    return Rational(QQ(x,preserve_ordering=preserve_ordering))
end



################################################################################
#
#  Copying
#
################################################################################

Base.copy(a::TropicalSemiringElem) = a

function Base.deepcopy_internal(x::TropicalSemiringElem, dict::IdDict)
  isinf(x) ? (return inf(parent(x))) : (return TropicalSemiringElem(x.parent, Base.deepcopy_internal(data(x), dict)))
end



################################################################################
#
#  Printing
#
################################################################################

# Hook into the fancy printing, we use (x) for finite values and ±∞ for infinity.
function AbstractAlgebra.expressify(x::TropicalSemiringElem{minOrMax}; context = nothing) where {minOrMax<:Union{typeof(min),typeof(max)}}
    if isinf(x)
        if Oscar.is_unicode_allowed()
            return minOrMax==typeof(min) ? "∞" : "-∞"
        else
            return minOrMax==typeof(min) ? "infty" : "-infty"
        end
    end
    return Expr(:call, "", expressify(data(x), context = context))
end

AbstractAlgebra.expressify(R::TropicalSemiring{typeof(min)}; context = nothing) = "Min tropical semiring"
AbstractAlgebra.expressify(R::TropicalSemiring{typeof(max)}; context = nothing) = "Max tropical semiring"
@enable_all_show_via_expressify TropicalSemiringElem
@enable_all_show_via_expressify TropicalSemiring



################################################################################
#
#  Equality and hash
#
################################################################################

function Base.:(==)(x::TropicalSemiringElem, y::TropicalSemiringElem)
  (isinf(x) && isinf(y)) && return true
  ((isinf(x) && !isinf(y)) || (!isinf(x) && isinf(y))) && return false
  return data(x) == data(y)
end

function Base.hash(x::TropicalSemiringElem, h::UInt)
  b = 0x4df38853cc07aa27   % UInt
  h = (isinf(x) ? hash(isinf(x), h) : hash(data(x), h))
  return xor(h, b)
end



################################################################################
#
#  Predicates
#  (see also isinf in basic access)
#
################################################################################

iszero(x::TropicalSemiringElem) = isinf(x)

isone(x::TropicalSemiringElem) = !isinf(x) && iszero(data(x))



################################################################################
#
#  Comparison
#  * min-plus semiring is ordered as usual: -inf < ... < -1 < 0 < 1 < ...
#  * max-plax semiring is ordered in reverse:       ... > -1 > 0 > 1 > ... > inf
#  (see Section 2.7 in Joswig: "Essentials of Tropical Combinatorics")
#
################################################################################

function isless(x::TropicalSemiringElem{typeof(min)}, y::TropicalSemiringElem{typeof(min)})
  iszero(x) && return false # x=-inf, no y is smaller
  iszero(y) && return true  # y=-inf, smaller than all x except x=-inf, which was handled above
  return data(x) < data(y)
end

function isless(x::TropicalSemiringElem{typeof(max)}, y::TropicalSemiringElem{typeof(max)})
  iszero(x) && return false # x=inf, no y is smaller
  iszero(y) && return true  # y=inf, smaller than all x except x=inf, which was handled above
  return data(x) > data(y)
end



################################################################################
#
#  Arithmetics
#
################################################################################

function Base.:(+)(x::TropicalSemiringElem{minOrMax}, y::TropicalSemiringElem{minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    iszero(x) && return deepcopy(y)                           # if x is zero, return y
    iszero(y) && return deepcopy(x)                           # if y is zero, return x
    return parent(x)(convention(parent(x))(data(x), data(y))) # otherwise, return their min / max
end

function Base.:(-)(x::TropicalSemiringElem, y::TropicalSemiringElem...)
  error("Tropical subtraction not defined (use tropical division for classical subtraction)")
end

function Base.:(*)(x::TropicalSemiringElem{minOrMax}, y::TropicalSemiringElem{minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    iszero(x) && return x # if x is zero, return it
    iszero(y) && return y # if y is zero, return it
    return parent(x)(data(x) + data(y)) # otherwise, return their sum
end

function divexact(x::TropicalSemiringElem{minOrMax}, y::TropicalSemiringElem{minOrMax}; check::Bool=true) where {minOrMax<:Union{typeof(min),typeof(max)}}
    @req !iszero(y) "dividing by (tropical) zero"
    return (iszero(x) ? x : parent(x)(data(x)-data(y)))
end

function Base.:(//)(x::TropicalSemiringElem{minOrMax}, y::TropicalSemiringElem{minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    return divexact(x,y)
end

function Base.:(/)(x::TropicalSemiringElem{minOrMax}, y::TropicalSemiringElem{minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    return divexact(x,y)
end

function inv(a::TropicalSemiringElem)
    @req !iszero(a) "inverting (tropical) zero"
    return parent(a)(-data(a))
end

function Base.:(^)(a::TropicalSemiringElem, n::QQFieldElem)
    if iszero(a)
        @req n>0 "dividing by (tropical) zero"
        return zero(parent(a))  # if a is zero, return zero
    end
    return parent(a)(data(a)*n) # otherwise (rational) multiply a by n
end

function Base.:(^)(a::TropicalSemiringElem, n::ZZRingElem)
    if iszero(a)
        @req n>0 "dividing by (tropical) zero"
        return zero(parent(a))  # if a is zero, return zero
    end
    return parent(a)(data(a)*n) # otherwise (rational) multiply a by n
end

function Base.:(^)(a::TropicalSemiringElem, n::Rational)
    if iszero(a)
        @req n>0 "dividing by (tropical) zero"
        return zero(parent(a))  # if a is zero, return zero
    end
    return parent(a)(data(a)*n) # otherwise (rational) multiply a by n
end

function Base.:(^)(a::TropicalSemiringElem, n::Integer)
    if iszero(a)
        @req n>0 "dividing by (tropical) zero"
        return zero(parent(a))  # if a is zero, return zero
    end
    return parent(a)(data(a)*n) # otherwise (rational) multiply a by n
end


################################################################################
#
#  Properties
#
################################################################################

function characteristic(::TropicalSemiring)
    error("characteristic of tropical semirings not supported")
end



################################################################################
#
#  helpers for polymake conversion
#
################################################################################

Polymake.convert_to_pm_type(::Type{<:TropicalSemiringElem{A}}) where A = Polymake.TropicalNumber{Polymake.convert_to_pm_type(A),Polymake.Rational}

function Base.convert(::Type{<:Polymake.TropicalNumber{PA}}, t::TropicalSemiringElem{A}) where {A <: Union{typeof(min),typeof(max)}, PA <: Union{Polymake.Min, Polymake.Max}}
  @req PA == Polymake.convert_to_pm_type(A) "cannot convert between different tropical conventions"
  isinf(t) ? Polymake.TropicalNumber{PA}() : Polymake.TropicalNumber{PA}(convert(Polymake.Rational, data(t)))
end

function (T::TropicalSemiring{A})(t::Polymake.TropicalNumber{PA}) where {A <: Union{typeof(min),typeof(max)}, PA <: Union{Polymake.Min, Polymake.Max}}
  @req PA == Polymake.convert_to_pm_type(A) "cannot convert between different tropical conventions"
  t == Polymake.zero(t) ? zero(T) : T(Polymake.scalar(t))
end
