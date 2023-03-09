export singular_coeff_ring, MPolyQuoRing, MPolyQuoRingElem, MPolyQuoIdeal
export quo, base_ring, modulus, gens, ngens, dim, simplify, default_ordering
export is_subset
export saturated_ideal
##############################################################################
#
# quotient rings
#
##############################################################################

@attributes mutable struct MPolyQuoRing{S} <: AbstractAlgebra.Ring
  I::MPolyIdeal{S}
  ordering::MonomialOrdering
  SQR::Singular.PolyRing # Singular quotient ring
  SQRGB::Singular.sideal # Singular Groebner basis defining quotient ring

  #= ordering, gb assure =#
  #= the current design decision is to fix the ordering to be default_ordering(R) =#
  #= that is, the user is not allowed to change this with the outer constructor quo =#
  function MPolyQuoRing(R::MPolyRing, I::MPolyIdeal, ordering::MonomialOrdering = default_ordering(R))
    @assert base_ring(I) === R
    r = new{elem_type(R)}()
    r.I = I
    r.ordering = ordering
    return r
  end
end

function groebner_assure(r::MPolyQuoRing)
  if isdefined(r, :SQRGB)
    return true
  end
  ordering = r.ordering
  groebner_assure(r.I, ordering)
  oscar_assure(r.I.gb[ordering])
  singular_assure(r.I.gb[ordering])
  SG = r.I.gb[ordering].gens.S
  r.SQR   = Singular.create_ring_from_singular_ring(Singular.libSingular.rQuotientRing(SG.ptr, base_ring(SG).ptr))
  r.SQRGB = Singular.Ideal(r.SQR, [r.SQR(0)])
  return true
end

function show(io::IO, Q::MPolyQuoRing)
  Hecke.@show_name(io, Q)
  Hecke.@show_special(io, Q)
  io = IOContext(io, :compact => true)
  print(io, "Quotient of $(base_ring(Q)) by $(modulus(Q))")
end

gens(Q::MPolyQuoRing) = [Q(x) for x = gens(base_ring(Q))]
ngens(Q::MPolyQuoRing) = ngens(base_ring(Q))
gen(Q::MPolyQuoRing, i::Int) = Q(gen(base_ring(Q), i))
Base.getindex(Q::MPolyQuoRing, i::Int) = Q(base_ring(Q)[i])::elem_type(Q)
base_ring(Q::MPolyQuoRing) = base_ring(Q.I)
coefficient_ring(Q::MPolyQuoRing) = coefficient_ring(base_ring(Q))
modulus(Q::MPolyQuoRing) = Q.I
oscar_groebner_basis(Q::MPolyQuoRing) = groebner_assure(Q) && return Q.I.gb[Q.ordering].O
singular_quotient_groebner_basis(Q::MPolyQuoRing) = groebner_assure(Q) && return Q.SQRGB
singular_origin_groebner_basis(Q::MPolyQuoRing) = groebner_assure(Q) && Q.I.gb[Q.ordering].gens.S
singular_quotient_ring(Q::MPolyQuoRing) = groebner_assure(Q) && Q.SQR
singular_poly_ring(Q::MPolyQuoRing) = singular_quotient_ring(Q)
singular_origin_ring(Q::MPolyQuoRing) = base_ring(singular_origin_groebner_basis(Q))
oscar_origin_ring(Q::MPolyQuoRing) = base_ring(Q)

default_ordering(Q::MPolyQuoRing) = default_ordering(base_ring(Q))

##############################################################################
#
# Quotient ring elements
#
##############################################################################

#TODO: think: do we want/ need to keep f on the Singular side to avoid conversions?
#      or use Bill's divrem to speed things up?
mutable struct MPolyQuoRingElem{S} <: RingElem
  f::S
  P::MPolyQuoRing{S}

  function MPolyQuoRingElem(f::S, P::MPolyQuoRing{S}) where {S}
    @assert parent(f) === base_ring(P)
    return new{S}(f, P)
  end
end

@enable_all_show_via_expressify MPolyQuoRingElem

function AbstractAlgebra.expressify(a::MPolyQuoRingElem; context = nothing)
  return expressify(a.f, context = context)
end

function Base.deepcopy_internal(a::MPolyQuoRingElem, dict::IdDict)
  return MPolyQuoRingElem(Base.deepcopy_internal(a.f, dict), a.P)
end

##############################################################################
#
# Quotient ring ideals
#
##############################################################################

# For ideals over quotient rings, we would like to delay the expensive
# construction of the singular quotient ring until the user does an operation
# that actually requires it.

@attributes mutable struct MPolyQuoIdeal{T} <: Ideal{T}
  gens::IdealGens{T}
  dim::Int
  gb::IdealGens{T}
  qRing::MPolyQuoRing

  function MPolyQuoIdeal(Ox::MPolyQuoRing{T}, si::Singular.sideal) where T <: MPolyRingElem
    singular_quotient_ring(Ox) == base_ring(si) || error("base rings must match")
    r = new{T}()
    r.gens = IdealGens(Ox, si)
    r.qRing = Ox
    br = base_ring(Ox)
    r.gens.gens.O = [br(g) for g = gens(r.gens.gens.S)]
    r.dim = -1
    return r
  end

  function MPolyQuoIdeal(Ox::MPolyQuoRing{T}, I::MPolyIdeal{T}) where T <: MPolyRingElem
    base_ring(Ox) === base_ring(I) || error("base rings must match")
    r = new{T}()
    r.gens = IdealGens(Ox, gens(I))
    r.qRing = Ox
    r.dim = -1
    return r
  end
  function MPolyQuoIdeal(Ox::MPolyQuoRing{T}, V::Vector{T}) where T <: MPolyRingElem
    r = new{T}()
    r.gens = IdealGens(Ox, V)
    r.qRing = Ox
    r.dim = -1
    return r
  end
end
@enable_all_show_via_expressify MPolyQuoIdeal

function AbstractAlgebra.expressify(a::MPolyQuoIdeal; context = nothing)
  return Expr(:call, :ideal, [expressify(g, context = context) for g in gens(a)]...)
end


@doc Markdown.doc"""
    base_ring(a::MPolyQuoIdeal)

Return the ambient ring of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> Q, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> a = ideal(Q, [x, y])
ideal(x, y)

julia> base_ring(a)
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z)
```
"""
function base_ring(a::MPolyQuoIdeal)
  return a.qRing
end

function oscar_assure(a::MPolyQuoIdeal)
  if isdefined(a.gens.gens, :O)
    return a.gens.gens.O
  end
  r = base_ring(base_ring(a))
  a.gens.gens.O = [r(g) for g = gens(a.gens.gens.S)]
end

function singular_assure(a::MPolyQuoIdeal)
  if isdefined(a.gens.gens, :S)
    return a.gens.S
  end
  a.gens.Sx = singular_poly_ring(base_ring(a))
  a.gens.S  = Singular.Ideal(a.gens.Sx, (a.gens.Sx).(gens(a)))
end

function groebner_assure(a::MPolyQuoIdeal)
  if !isdefined(a, :gb)
    singular_assure(a)
    a.gb = IdealGens(base_ring(a), Singular.std(a.gens.S))
    a.gb.gens.S.isGB = a.gb.isGB = true
  end
end


@doc Markdown.doc"""
    gens(a::MPolyQuoIdeal)

Return the generators of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, QQMPolyRingElem[x, y, z])

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]))
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z) defined by a julia-function with inverse)

julia> a = ideal(A, [x-y])
ideal(x - y)

julia> gens(a)
1-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x - y
```
"""
function gens(a::MPolyQuoIdeal)
  oscar_assure(a)
  return map(a.gens.Ox, a.gens.O)
end

gen(a::MPolyQuoIdeal, i::Int) = gens(a)[i]
getindex(a::MPolyQuoIdeal, i::Int) = gen(a, i)

@doc Markdown.doc"""
    ngens(a::MPolyQuoIdeal)

Return the number of generators of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, QQMPolyRingElem[x, y, z])

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]))
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z) defined by a julia-function with inverse)

julia> a = ideal(A, [x-y])
ideal(x - y)

julia> ngens(a)
1
```
"""
function ngens(a::MPolyQuoIdeal)
  return length(gens(a))
end


# powers, addition and multiplication do not require the singular quotient ring

@doc Markdown.doc"""
    :^(a::MPolyQuoIdeal, m::Int)

Return the `m`-th power of `a`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, [x^2-y, y^2-x+y]);

julia> a = ideal(A, [x+y])
ideal(x + y)

julia> a^2
ideal(x^2 + 2*x*y + y^2)
```
"""
function Base.:^(a::MPolyQuoIdeal, m::Int)
  singular_assure(a)
  return MPolyQuoIdeal(base_ring(a), a.gens.S^m)
end

@doc Markdown.doc"""
    :+(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return the sum of `a` and `b`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, [x^2-y, y^2-x+y]);

julia> a = ideal(A, [x+y])
ideal(x + y)

julia> b = ideal(A, [x^2+y^2, x+y])
ideal(x^2 + y^2, x + y)

julia> a+b
ideal(x + y, x^2 + y^2)
```
"""
function Base.:+(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  base_ring(a) == base_ring(b) || error("base rings must match")
  singular_assure(a)
  singular_assure(b)
  return MPolyQuoIdeal(base_ring(a), a.gens.S + b.gens.S)
end

@doc Markdown.doc"""
    :*(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return the product of `a` and `b`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, [x^2-y, y^2-x+y]);

julia> a = ideal(A, [x+y])
ideal(x + y)

julia> b = ideal(A, [x^2+y^2, x+y])
ideal(x^2 + y^2, x + y)

julia> a*b
ideal(x^3 + x^2*y + x*y^2 + y^3, x^2 + 2*x*y + y^2)
```
"""
function Base.:*(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  base_ring(a) == base_ring(b) || error("base rings must match")
  singular_assure(a)
  singular_assure(b)
  return MPolyQuoIdeal(base_ring(a), a.gens.S * b.gens.S)
end

@doc Markdown.doc"""
    intersect(a::MPolyQuoIdeal{T}, bs::MPolyQuoIdeal{T}...) where T

Return the intersection of two or more ideals.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [x^2-y^3, x-y]));

julia> a = ideal(A, [y^2])
ideal(y^2)

julia> b = ideal(A, [x])
ideal(x)

julia> intersect(a,b)
ideal(x*y)
```
"""
function intersect(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}...) where T
  singular_assure(a)
  as = a.gens.S
  for g in b
    base_ring(g) == base_ring(a) || error("base rings must match")
    singular_assure(g)
    gs = g.gens.S
    as = Singular.intersection(as, gs)
  end
  return MPolyQuoIdeal(base_ring(a), as)
end

#######################################################

@doc Markdown.doc"""
    quotient(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return the ideal quotient of `a` by `b`. Alternatively, use `a:b`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [x^2-y^3, x-y]));

julia> a = ideal(A, [y^2])
ideal(y^2)

julia> b = ideal(A, [x])
ideal(x)

julia> a:b
ideal(y)
```
"""
function quotient(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  base_ring(a) == base_ring(b) || error("base rings must match")

  singular_assure(a)
  singular_assure(b)
  return MPolyQuoIdeal(base_ring(a), Singular.quotient(a.gens.S, b.gens.S))
end
(::Colon)(a::MPolyQuoIdeal, b::MPolyQuoIdeal) = quotient(a, b)

# TODO: replace by a more efficient method!
@attr function is_prime(I::MPolyQuoIdeal)
  return is_prime(saturated_ideal(I))
end

# The following is to streamline the programmer's
# interface for the use of the four standard rings
# for the schemes `MPolyRing`, `MPolyQuoRing`, `MPolyLocRing`,
# and `MPolyQuoLocRing` together with their ideals.
# We return the preimage of the given ideal under the
# canonical map from the underlying free polynomial ring.
@attr function saturated_ideal(I::MPolyQuoIdeal)
  R = base_ring(base_ring(I))
  J = ideal(R, lift.(gens(I))) + modulus(base_ring(I))
  return J
end

@doc Markdown.doc"""
    is_zero(a::MPolyQuoIdeal)

Return `true` if `a` is the zero ideal, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, [x^2-y, y^2-x+y]);

julia> a = ideal(A, [x^2+y^2, x+y])
ideal(x^2 + y^2, x + y)

julia> is_zero(a)
false

julia> b = ideal(A, [x^2-y])
ideal(x^2 - y)

julia> is_zero(b)
true
```
"""
function is_zero(a::MPolyQuoIdeal)
  R = base_ring(a)
  singular_assure(a)
  return Singular.iszero(Singular.reduce(a.gens.S, singular_quotient_groebner_basis(R)))
end

@doc Markdown.doc"""
    ideal(A::MPolyQuoRing{T}, V::Vector{T}) where T <: MPolyRingElem

Given a (graded) quotient ring `A=R/I` and a vector `V` of (homogeneous) polynomials in `R`,
create the ideal of `A` which is generated by the images of the entries of `V`.

    ideal(A::MPolyQuoRing{T}, V::Vector{MPolyQuoRingElem{T}}) where T <: MPolyRingElem

Given a (graded) quotient ring `A` and a vector `V` of (homogeneous) elements of `A`,
create the ideal of `A` which is generated by the entries of `V`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [x^2-y^3, x-y]));

julia> I = ideal(A, [x^2-y])
ideal(x^2 - y)
```
```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> B, _ = quo(S, ideal(S, [x^2*z-y^3, x-y]));

julia> J = ideal(B, [x^2-y^2])
ideal(x^2 - y^2)
```
"""
function ideal(A::MPolyQuoRing{T}, V::Vector{T}) where T <: MPolyRingElem
  #@assert length(V) > 0
  if length(V) == 0
    return MPolyQuoIdeal(A, elem_type(base_ring(A))[])
  end
  for p in V
    base_ring(A) == parent(p) || error("parents must match")
  end
  return MPolyQuoIdeal(A, V)
end
function ideal(A::MPolyQuoRing{T}, V::Vector{MPolyQuoRingElem{T}}) where T <: MPolyRingElem
  #@assert length(V) > 0
  if length(V) == 0
    return MPolyQuoIdeal(A, ideal(base_ring(A), elem_type(base_ring(A))[]))
  end
  for p in V
    A == parent(p) || error("parents must match")
  end
  return MPolyQuoIdeal(A, ideal(base_ring(A), map(p->p.f, V)))
end

function ideal(A::MPolyQuoRing{T}, x::T) where T <: MPolyRingElem
  return ideal(A,[x])
end

function ideal(A::MPolyQuoRing{T}, x::MPolyQuoRingElem{T}) where T <: MPolyRingElem
  return ideal(A,[x])
end
##################################################################

parent_type(::MPolyQuoRingElem{S}) where S = MPolyQuoRing{S}
parent_type(::Type{MPolyQuoRingElem{S}}) where S = MPolyQuoRing{S}
elem_type(::MPolyQuoRing{S})  where S= MPolyQuoRingElem{S}
elem_type(::Type{MPolyQuoRing{S}})  where S= MPolyQuoRingElem{S}

canonical_unit(a::MPolyQuoRingElem) = one(parent(a))

parent(a::MPolyQuoRingElem) = a.P

function check_parent(a::MPolyQuoRingElem, b::MPolyQuoRingElem)
  a.P == b.P || error("wrong parents")
  return true
end

+(a::MPolyQuoRingElem, b::MPolyQuoRingElem) = check_parent(a, b) && MPolyQuoRingElem(a.f+b.f, a.P)

-(a::MPolyQuoRingElem, b::MPolyQuoRingElem) = check_parent(a, b) && MPolyQuoRingElem(a.f-b.f, a.P)

-(a::MPolyQuoRingElem) = MPolyQuoRingElem(-a.f, a.P)

*(a::MPolyQuoRingElem, b::MPolyQuoRingElem) = check_parent(a, b) && simplify(MPolyQuoRingElem(a.f*b.f, a.P))

^(a::MPolyQuoRingElem, b::Base.Integer) = simplify(MPolyQuoRingElem(Base.power_by_squaring(a.f, b), a.P))

*(a::MPolyQuoRingElem, b::QQFieldElem) = simplify(MPolyQuoRingElem(a.f * b, a.P))

*(a::MPolyQuoRingElem, b::ZZRingElem) = simplify(MPolyQuoRingElem(a.f * b, a.P))

*(a::QQFieldElem, b::MPolyQuoRingElem) = simplify(MPolyQuoRingElem(a * b.f, b.P))

*(a::ZZRingElem, b::MPolyQuoRingElem) = simplify(MPolyQuoRingElem(a * b.f, b.P))

#*(a::MPolyQuoRingElem, b::MPolyQuoRingElem) = check_parent(a, b) && MPolyQuoRingElem(a.f*b.f, a.P)
#
#^(a::MPolyQuoRingElem, b::Base.Integer) = MPolyQuoRingElem(Base.power_by_squaring(a.f, b), a.P)

function Oscar.mul!(a::MPolyQuoRingElem, b::MPolyQuoRingElem, c::MPolyQuoRingElem)
  a.f = b.f*c.f
  return a
end

function Oscar.addeq!(a::MPolyQuoRingElem, b::MPolyQuoRingElem)
  a.f += b.f
  return a
end

@doc Markdown.doc"""
    simplify(a::MPolyQuoIdeal)

If `a` is an ideal of the quotient of a multivariate polynomial ring `R` by an ideal `I` of `R`, say,
replace the internal polynomial representative of each generator of `a` by its normal form 
mod `I` with respect to the `default_ordering` on `R`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> a = ideal(A, [x^3*y^4-x+y, x*y+y^2*x])
ideal(x^3*y^4 - x + y, x*y^2 + x*y)

julia> gens(a)
2-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x^3*y^4 - x + y
 x*y^2 + x*y

julia> simplify(a)
ideal(x^2*y^3 - x + y, x*y^2 + x*y)

julia> gens(a)
2-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x^2*y^3 - x + y
 x*y^2 + x*y
```
"""
function simplify(a::MPolyQuoIdeal)
  Q = base_ring(a)
  R = base_ring(Q)
  singular_assure(a)
  red  = reduce(a.gens.S, singular_quotient_groebner_basis(Q))
  SQ   = singular_poly_ring(Q)
  si   = Singular.Ideal(SQ, unique!(gens(red)))
  a.gens.S = si
  a.gens.O = [R(g) for g = gens(a.gens.S)]

  return a
end

#######################################################

@doc Markdown.doc"""
    ideal_membership(f::MPolyQuoRingElem{T}, a::MPolyQuoIdeal{T}) where T

Return `true` if `f` is contained in `a`, `false` otherwise. Alternatively, use `f in a`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> a = ideal(A, [x^3*y^4-x+y, x*y+y^2*x])
ideal(x^3*y^4 - x + y, x*y^2 + x*y)

julia> f = A(x^2*y^3-x+y)
x^2*y^3 - x + y

julia> f in a
true
```
"""
function ideal_membership(a::MPolyQuoRingElem{T}, b::MPolyQuoIdeal{T}) where T
  parent(a) == base_ring(b) || error("base rings must match")
  groebner_assure(b)
  SR = singular_poly_ring(base_ring(b))
  as = simplify(a)
  return Singular.iszero(Singular.reduce(SR(as), b.gb.gens.S))
end

Base.:in(a::MPolyQuoRingElem, b::MPolyQuoIdeal) = ideal_membership(a, b)


@doc Markdown.doc"""
    is_subset(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return `true` if `a` is contained in `b`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> a = ideal(A, [x^3*y^4-x+y, x*y+y^2*x])
ideal(x^3*y^4 - x + y, x*y^2 + x*y)

julia> b = ideal(A, [x^3*y^3-x+y, x^2*y+y^2*x])
ideal(x^3*y^3 - x + y, x^2*y + x*y^2)

julia> is_subset(a,b)
false

julia> is_subset(b,a)
true
```
"""
function is_subset(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  base_ring(a) == base_ring(b) || error("base rings must match")
  as = simplify(a)
  groebner_assure(b)
  return Singular.iszero(Singular.reduce(as.gens.S, b.gb.gens.S))
end

@doc Markdown.doc"""
    ==(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return `true` if `a` is equal to `b`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> a = ideal(A, [x^3*y^4-x+y, x*y+y^2*x])
ideal(x^3*y^4 - x + y, x*y^2 + x*y)

julia> b = ideal(A, [x^3*y^3-x+y, x^2*y+y^2*x])
ideal(x^3*y^3 - x + y, x^2*y + x*y^2)

julia> a == b
false
```
"""
function Base.:(==)(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  return issubset(a, b) && issubset(b, a)
end

@doc Markdown.doc"""
    simplify(f::MPolyQuoRingElem)

If `f` is an element of the quotient of a multivariate polynomial ring `R` by an ideal `I` of `R`, say,
replace the internal polynomial representative of `f` by its normal form mod `I` with respect to 
the `default_ordering` on `R`.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, ["x"]);

julia> A, p = quo(R, ideal(R, [x^4]));

julia> f = p(2*x^6 + x^3 + x)
2*x^6 + x^3 + x

julia> simplify(f)
x^3 + x

julia> f
x^3 + x
```
"""
function simplify(f::MPolyQuoRingElem)
  R  = parent(f)
  OR = oscar_origin_ring(R)
  SR = singular_origin_ring(R)
  G  = singular_origin_groebner_basis(R)
  g  = f.f
  f.f = OR(reduce(SR(g), G))
  return f::elem_type(R)
end


@doc Markdown.doc"""
    ==(f::MPolyQuoRingElem{T}, g::MPolyQuoRingElem{T}) where T

Return `true` if `f` is equal to `g`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, ["x"]);

julia> A, p = quo(R, ideal(R, [x^4]));

julia> f = p(x-x^6)
-x^6 + x

julia> g = p(x)
x

julia> f == g
true
```
"""
function ==(f::MPolyQuoRingElem{T}, g::MPolyQuoRingElem{T}) where T
  check_parent(f, g)
  simplify(f)
  simplify(g)
  return f.f == g.f
end

@doc Markdown.doc"""
    quo(R::MPolyRing, I::MPolyIdeal) -> MPolyQuoRing, Map

Create the quotient ring `R/I` and return the new
ring as well as the projection map `R` $\to$ `R/I`.

    quo(R::MPolyRing, V::Vector{MPolyRingElem}) -> MPolyQuoRing, Map

As above, where `I` is the ideal of `R` generated by the polynomials in `V`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, p = quo(R, ideal(R, [x^2-y^3, x-y]));

julia> A
Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - y^3, x - y)

julia> typeof(A)
MPolyQuoRing{QQMPolyRingElem}

julia> typeof(x)
QQMPolyRingElem

julia> p
Map from
Multivariate Polynomial Ring in x, y over Rational Field to Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - y^3, x - y) defined by a julia-function with inverse

julia> p(x)
x

julia> typeof(p(x))
MPolyQuoRingElem{QQMPolyRingElem}
```
```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> B, _ = quo(S, ideal(S, [x^2*z-y^3, x-y]))
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] by ideal(x^2*z - y^3, x - y), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] by ideal(x^2*z - y^3, x - y) defined by a julia-function with inverse)

julia> typeof(B)
MPolyQuoRing{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}
```
"""
function quo(R::MPolyRing, I::MPolyIdeal)
  q = MPolyQuoRing(R, I)
  function im(a::MPolyRingElem)
    parent(a) !== R && error("Element not in the domain of the map")
    return MPolyQuoRingElem(a, q)
  end
  function pr(a::MPolyQuoRingElem)
    return a.f
  end
  return q, MapFromFunc(im, pr, R, q)
end

function quo(R::MPolyRing, I::Vector{<:MPolyRingElem})
  return quo(R, ideal(I))
end

function quo(R::MPolyRing, f::MPolyRingElem...)
  return quo(R, ideal(collect(f)))
end

lift(a::MPolyQuoRingElem) = a.f

(Q::MPolyQuoRing)() = MPolyQuoRingElem(base_ring(Q)(), Q)

function (Q::MPolyQuoRing)(a::MPolyQuoRingElem)
  parent(a) !== Q && error("Parent mismatch")
  return a
end

function (Q::MPolyQuoRing{S})(a::S) where {S <: MPolyRingElem}
  base_ring(Q) === parent(a) || error("Parent mismatch")
  return MPolyQuoRingElem(a, Q)
end

function (Q::MPolyQuoRing)(a::MPolyRingElem)
  return Q(base_ring(Q)(a))
end

function (Q::MPolyQuoRing)(a::Singular.spoly)
  @assert singular_poly_ring(Q) == parent(a)
  return MPolyQuoRingElem(base_ring(Q)(a), Q)
end

function (S::Singular.PolyRing)(a::MPolyQuoRingElem)
  Q = parent(a)
  @assert singular_poly_ring(Q) == S
  return S(a.f)
end

(Q::MPolyQuoRing)(a) = MPolyQuoRingElem(base_ring(Q)(a), Q)

zero(Q::MPolyQuoRing) = Q(0)
one(Q::MPolyQuoRing) = Q(1)

function is_invertible_with_inverse(a::MPolyQuoRingElem)
  # TODO:
  # Eventually, the code below should be replaced
  # by a call to `coordinates` over the ring `parent(a)`.
  # This should then use relative groebner bases and
  # make use of the caching of previously computed GBs
  # of the modulus of `parent(a)`.

  Q = parent(a)
  J = oscar_groebner_basis(Q)
  J = vcat(J, [a.f])
  j, T = standard_basis_with_transformation_matrix(ideal(J))
  if is_constant(j[1]) && is_unit(first(coefficients(j[1])))
    @assert ncols(T) == 1
    return true, inv(first(coefficients(j[1])))*Q(T[end, 1])
  end
  return false, a
end

is_unit(a::MPolyQuoRingElem) = is_invertible_with_inverse(a)[1]

function inv(a::MPolyQuoRingElem)
  fl, b = is_invertible_with_inverse(a)
  fl || error("Element not invertible")
  return b
end


"""
Converts a sparse-Singular vector of polynomials to an Oscar sparse row.
"""
function sparse_row(R::MPolyRing, M::Singular.svector{<:Singular.spoly})
  v = Dict{Int, MPolyBuildCtx}()
  for (i, e, c) = M
    if !haskey(v, i)
      v[i] = MPolyBuildCtx(R)
    end
    push_term!(v[i], base_ring(R)(c), e)
  end
  pos_value_vector::Vector{Tuple{Int, elem_type(R)}} = [(k,finish(v)) for (k,v) = v]
  return sparse_row(R, pos_value_vector)
end

"""
Converts a sparse-Singular vector of polynomials to an Oscar sparse row.
Collect only the column indices in `U`.
"""
function sparse_row(R::MPolyRing, M::Singular.svector{<:Singular.spoly}, U::UnitRange)
  v = Dict{Int, MPolyBuildCtx}()
  for (i, e, c) = M
    (i in U) || continue
    if !haskey(v, i)
      v[i] = MPolyBuildCtx(R)
    end
    push_term!(v[i], base_ring(R)(c), e)
  end
  pos_value_vector::Vector{Tuple{Int, elem_type(R)}} = [(k,finish(v)) for (k,v) = v]
  return sparse_row(R, pos_value_vector)
end

"""
Converts the sparse-Singular matrix (`Module`) row by row to an Oscar sparse-matrix.
Only the row indices (generators) in `V` and the column indices in `U` are converted.
"""
function sparse_matrix(R::MPolyRing, M::Singular.Module, V::UnitRange, U::UnitRange)
  S = sparse_matrix(R)
  for g = 1:Singular.ngens(M)
    (g in V) || continue
    push!(S, sparse_row(R, M[g], U))
  end
  return S
end

"""
Converts the sparse-Singular matrix (`Module`) row by row to an Oscar sparse-matrix.
"""
function sparse_matrix(R::MPolyRing, M::Singular.Module)
  S = sparse_matrix(R)
  for g = 1:Singular.ngens(M)
    push!(S, sparse_row(R, M[g]))
  end
  S.r = ngens(M)
  S.c = rank(M)
  return S
end

"""
Converts the sparse-Singular matrix (`Module`) row by row to an Oscar dense-matrix.
"""
function matrix(R::MPolyRing, M::Singular.Module)
  return matrix(sparse_matrix(R, M))
end

function divides(a::MPolyQuoRingElem, b::MPolyQuoRingElem)
  check_parent(a, b)
  iszero(a) && iszero(b) && return (true, zero(parent(a)))
  iszero(b) && error("cannot divide by zero")

  Q = parent(a)
  J = oscar_groebner_basis(Q)

  BS = IdealGens([a.f], keep_ordering = false)
  singular_assure(BS)

  J = vcat(J, [b.f])
  BJ = IdealGens(J, keep_ordering = false)
  singular_assure(BJ)

  s, rest = Singular.lift(BJ.S, BS.S)
  if !iszero(rest)
    return false, a
  end
  return true, Q(sparse_matrix(base_ring(Q), s, 1:1, length(J):length(J))[1, length(J)])
end

#TODO: find a more descriptive, meaningful name
function _kbase(Q::MPolyQuoRing)
  G = singular_origin_groebner_basis(Q)
  s = Singular.kbase(G)
  if iszero(s)
    error("ideal was no zero-dimensional")
  end
  return [base_ring(Q)(x) for x = gens(s)]
end

#TODO: the reverse map...
# problem: the "canonical" reps are not the monomials.
function vector_space(K::AbstractAlgebra.Field, Q::MPolyQuoRing)
  R = base_ring(Q)
  @assert K == base_ring(R)
  l = _kbase(Q)
  V = free_module(K, length(l))
  function im(a::Generic.FreeModuleElem)
    @assert parent(a) == V
    b = R(0)
    for i=1:length(l)
      c = a[i]
      if !iszero(c)
        b += c*l[i]
      end
    end
    return Q(b)
  end
  return V, MapFromFunc(im, V, Q)
end

# To fix printing of fraction fields of MPolyQuoRing
function AbstractAlgebra.expressify(a::AbstractAlgebra.Generic.Frac{T};
    context = nothing) where {T <: MPolyQuoRingElem}
  n = numerator(a, false)
  d = denominator(a, false)
  if isone(d)
    return expressify(n, context = context)
  else
    return Expr(:call, ://, expressify(n, context = context),
                expressify(d, context = context))
  end
end

################################################################################
#
#  Graded functionality
#
################################################################################

function grading(R::MPolyQuoRing)
  if base_ring(R) isa MPolyDecRing
    return grading(base_ring(R))
  else
    error("Underlying polynomial ring must be graded")
  end
end

@doc Markdown.doc"""
    degree(f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given a homogeneous element `f` of a graded affine algebra, return the degree of `f`.

    degree(::Type{Vector{Int}}, f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given a homogeneous element `f` of a $\mathbb Z^m$-graded affine algebra, return the degree of `f`, converted to a vector of integer numbers.

    degree(::Type{Int}, f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given a homogeneous element `f` of a $\mathbb Z$-graded affine algebra, return the degree of `f`, converted to an integer number.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"] );

julia> A, p = quo(R, ideal(R, [y-x, z^3-x^3]))
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] by ideal(-x + y, -x^3 + z^3), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] by ideal(-x + y, -x^3 + z^3) defined by a julia-function with inverse)

julia> f = p(y^2-x^2+z^4)
-x^2 + y^2 + z^4

julia> degree(f)
graded by [4]

julia> typeof(degree(f))
GrpAbFinGenElem

julia> degree(Int, f)
4

julia> typeof(degree(Int, f))
Int64
```
"""
function degree(a::MPolyQuoRingElem{<:MPolyDecRingElem})
  simplify(a)
  @req !iszero(a) "Element must be non-zero"
  return degree(a.f)
end

function degree(::Type{Int}, a::MPolyQuoRingElem{<:MPolyDecRingElem})
  @assert is_z_graded(base_ring(parent(a)))
  return Int(degree(a)[1])
end

function degree(::Type{Vector{Int}}, a::MPolyQuoRingElem{<:MPolyDecRingElem})
  @assert is_zm_graded((base_ring(parent(a))))
  d = degree(a)
  return Int[d[i] for i=1:ngens(parent(d))]
end

is_filtered(q::MPolyQuoRing) = is_filtered(base_ring(q))
is_graded(q::MPolyQuoRing) = is_graded(base_ring(q))

@doc Markdown.doc"""
    homogeneous_component(f::MPolyQuoRingElem{<:MPolyDecRingElem}, g::GrpAbFinGenElem)

Given an element `f` of a graded affine algebra, and given an element `g` of the
grading group of that algebra, return the homogeneous component of `f` of degree `g`.

    homogeneous_component(f::MPolyQuoRingElem{<:MPolyDecRingElem}, g::Vector{<:IntegerUnion})

Given an element `f` of a $\mathbb  Z^m$-graded affine algebra `A`, say, and given
a vector `g` of $m$ integers, convert `g` into an element of the grading group of `A`,
and return the homogeneous component of `f` whose degree is that element.

    homogeneous_component(f::MPolyQuoRingElem{<:MPolyDecRingElem}, g::IntegerUnion)

Given an element `f` of a $\mathbb  Z$-graded affine algebra `A`, say, and given
an integer `g`, convert `g` into an element of the grading group of `A`,
and return the homogeneous component of `f` whose degree is that element.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, p = quo(R, ideal(R, [y-x, z^3-x^3]));

julia> f = p(y^2-x^2+x*y*z+z^4)
-x^2 + x*y*z + y^2 + z^4

julia> homogeneous_component(f, 4)
z^4
```
"""
function homogeneous_component(a::MPolyQuoRingElem{<:MPolyDecRingElem}, d::GrpAbFinGenElem)
  simplify(a)
  return homogeneous_component(a.f, d)
end

function homogeneous_component(a::MPolyQuoRingElem{<:MPolyDecRingElem}, g::IntegerUnion)
  @assert is_z_graded(base_ring(parent(a)))
  return homogeneous_component(a, grading_group(base_ring(parent(a)))([g]))
end

function homogeneous_component(a::MPolyQuoRingElem{<:MPolyDecRingElem}, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(base_ring(parent(a)))
  return homogeneous_component(a, grading_group(base_ring(parent(a)))(g))
end

@doc Markdown.doc"""
    homogeneous_components(f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given an element `f` of a graded affine algebra, return the homogeneous components of `f`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, p = quo(R, ideal(R, [y-x, z^3-x^3]));

julia> f = p(y^2-x^2+x*y*z+z^4)
-x^2 + x*y*z + y^2 + z^4

julia> homogeneous_components(f)
Dict{GrpAbFinGenElem, MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}} with 2 entries:
  [4] => z^4
  [3] => y^2*z
```
"""
function homogeneous_components(a::MPolyQuoRingElem{<:MPolyDecRingElem})
  simplify(a)
  h = homogeneous_components(a.f)
  return Dict{keytype(h), typeof(a)}(x => parent(a)(y) for (x, y) in h)
end

@doc Markdown.doc"""
    is_homogeneous(f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given an element `f` of a graded affine algebra, return `true` if `f` is homogeneous, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, p = quo(R, ideal(R, [y-x, z^3-x^3]));

julia> f = p(y^2-x^2+z^4)
-x^2 + y^2 + z^4

julia> is_homogeneous(f)
true

julia> f
z^4
```
"""
function is_homogeneous(a::MPolyQuoRingElem{<:MPolyDecRingElem})
  simplify(a)
  return is_homogeneous(a.f)
end

@doc Markdown.doc"""
    grading_group(A::MPolyQuoRing{<:MPolyDecRingElem})

If `A` is, say, `G`-graded, return `G`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [x^2*z-y^3, x-y]));

julia> grading_group(A)
GrpAb: Z
```
"""
function grading_group(A::MPolyQuoRing{<:MPolyDecRingElem})
  return grading_group(base_ring(A))
end

function hash(w::MPolyQuoRingElem, u::UInt)
  simplify(w)
  return hash(w.f, u)
end

function homogeneous_component(W::MPolyQuoRing{<:MPolyDecRingElem}, d::GrpAbFinGenElem)
  #TODO: lazy: ie. no enumeration of points
  #      apparently it is possible to get the number of points faster than the points
  D = parent(d)
  @assert D == grading_group(W)
  R = base_ring(W)

  H, mH = homogeneous_component(R, d)
  B = Set{elem_type(W)}()
  for h = basis(H)
    b = W(mH(h))
    if !iszero(b)
      push!(B, b)
    end
  end
  B = [x for x = B]

  M, h = vector_space(base_ring(R), B, target = W)
  set_attribute!(M, :show => show_homo_comp, :data => (W, d))
  return M, h
end

@doc Markdown.doc"""
    dim(a::MPolyQuoIdeal)

Return the Krull dimension of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> a = ideal(A, [x-y])
ideal(x - y)

julia> dim(a)
0
```
"""
function dim(a::MPolyQuoIdeal)
  if a.dim > -1
    return a.dim
  end
  groebner_assure(a)
  a.dim = Singular.dimension(a.gb.S)
  return a.dim
end

##################################
### Tests on graded quotient rings
##################################

function is_standard_graded(A::MPolyQuoRing)
  return is_standard_graded(base_ring(A))
end

function is_z_graded(A::MPolyQuoRing)
  return is_z_graded(base_ring(A))
end

function is_zm_graded(A::MPolyQuoRing)
  return is_zm_graded(base_ring(A))
end

function is_positively_graded(A::MPolyQuoRing)
  return is_positively_graded(base_ring(A))
end

##################################
#######################################################
@doc Markdown.doc"""
    minimal_generating_set(I::MPolyQuoIdeal{<:MPolyDecRingElem})

Given a homogeneous ideal `a` of a graded affine algebra over a field,
return an array containing a minimal set of generators of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> V = [x, z^2, x^3+y^3, y^4, y*z^5];

julia> I = ideal(R, V)
ideal(x, z^2, x^3 + y^3, y^4, y*z^5)

julia> A, p = quo(R, ideal(R, [x-y]));

julia> a = ideal(A, [p(x) for x in V]);

julia> minimal_generating_set(a)
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 x
 z^2
```
"""
function minimal_generating_set(I::MPolyQuoIdeal{<:MPolyDecRingElem}; ordering::MonomialOrdering = default_ordering(base_ring(base_ring(I))))
  # This only works / makes sense for homogeneous ideals. So far ideals in an
  # MPolyDecRing are forced to be homogeneous though.

  Q = base_ring(I)

  @assert is_graded(Q)

  if !(coefficient_ring(Q) isa AbstractAlgebra.Field)
    throw(ArgumentError("The coefficient ring must be a field."))
  end

  QS = singular_poly_ring(Q)
  singular_assure(I)

  IS = I.gens.S
  GC.@preserve IS QS begin
    ptr = Singular.libSingular.idMinBase(IS.ptr, QS.ptr)
    gensS = gens(typeof(IS)(QS, ptr))
  end

  i = 1
  while i <= length(gensS)
    if iszero(gensS[i])
      deleteat!(gensS, i)
    else
      i += 1
    end
  end

  return elem_type(Q)[ Q(f) for f in gensS ]
end

################################################################################
#
#  Promote rule
#
################################################################################

function AbstractAlgebra.promote_rule(::Type{MPolyQuoRingElem{S}}, ::Type{T}) where {S, T <: RingElem}
  if AbstractAlgebra.promote_rule(S, T) === S
    return MPolyQuoRingElem{S}
  else
    return Union{}
  end
end

@attr function _is_integral_domain(A::MPolyQuoRing)
  return is_prime(modulus(A))
end

