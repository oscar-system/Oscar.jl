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

  function MPolyQuoRing(R::MPolyRing, I::MPolyIdeal, ordering::MonomialOrdering = default_ordering(R))
    @req base_ring(I) === R "Base rings must be the same."
    @req is_global(ordering) "Ordering must be global."
    r = new{elem_type(R)}()
    r.I = I
    r.ordering = ordering
    return r
  end
end

function _groebner_basis(r::MPolyQuoRing)
  # Assure that the correspondence to Singular is set up
  # Note that this is the only way to create the Singular quotient ring:
  # 1. Translate the modulus to a Singular ideal
  # 2. Create the quotient ring (rQuotientRing) on the Singular side
  # 3. Create a Singular.PolyRing (create_ring_from_singular_ring).
  if isdefined(r, :SQRGB)
    return true
  end
  ordering = r.ordering
  groebner_basis(r.I, ordering=ordering)
  oscar_assure(r.I.gb[ordering])
  SG = singular_generators(r.I.gb[ordering], ordering)
  r.SQR   = Singular.create_ring_from_singular_ring(Singular.libSingular.rQuotientRing(SG.ptr, base_ring(SG).ptr))
  r.SQRGB = Singular.Ideal(r.SQR, [r.SQR(0)])
  return true
end

function Base.show(io::IO, ::MIME"text/plain", Q::MPolyQuoRing)
   io = pretty(io)
   println(io, "Quotient")
   print(io, Indent(), "of ", Lowercase())
   show(io, MIME("text/plain"), base_ring(Q))
   println(io)
   print(io, "by ", Lowercase(), modulus(Q))
   print(io, Dedent())
end

function Base.show(io::IO, Q::MPolyQuoRing)
  @show_name(io, Q)
  @show_special(io, Q)
  if is_terse(io)
    # no nested printing
    print(io, "Quotient of multivariate polynomial ring")
  else
    # nested printing allowed
    io = pretty(io)
    print(io, "Quotient of multivariate polynomial ring by ", Lowercase(), modulus(Q))
  end
end

gens(Q::MPolyQuoRing) = [Q(x) for x = gens(base_ring(Q))]
number_of_generators(Q::MPolyQuoRing) = number_of_generators(base_ring(Q))
gen(Q::MPolyQuoRing, i::Int) = Q(gen(base_ring(Q), i))
base_ring(Q::MPolyQuoRing) = base_ring(Q.I)
base_ring_type(::Type{MPolyQuoRing{S}}) where {S} = base_ring_type(S)
coefficient_ring(Q::MPolyQuoRing) = coefficient_ring(base_ring(Q))
modulus(Q::MPolyQuoRing) = Q.I
oscar_groebner_basis(Q::MPolyQuoRing) = _groebner_basis(Q) && return Q.I.gb[Q.ordering].O
singular_quotient_groebner_basis(Q::MPolyQuoRing) = _groebner_basis(Q) && return Q.SQRGB
singular_origin_groebner_basis(Q::MPolyQuoRing) = _groebner_basis(Q) && Q.I.gb[Q.ordering].gens.S
singular_quotient_ring(Q::MPolyQuoRing) = _groebner_basis(Q) && Q.SQR
singular_poly_ring(Q::MPolyQuoRing; keep_ordering::Bool = false) = singular_quotient_ring(Q)
singular_poly_ring(Q::MPolyQuoRing, ordering::MonomialOrdering) = singular_quotient_ring(Q)
singular_origin_ring(Q::MPolyQuoRing) = base_ring(singular_origin_groebner_basis(Q))
oscar_origin_ring(Q::MPolyQuoRing) = base_ring(Q)

default_ordering(Q::MPolyQuoRing) = default_ordering(base_ring(Q))

# Only for fields for now because of things like char(ZZ[x, y]/<2>) = 2
function characteristic(Q::MPolyQuoRing{<:MPolyRingElem{T}}) where {T <: FieldElement}
  if is_zero(one(Q))
    return 1
  end
  return characteristic(coefficient_ring(Q))
end

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
  simplified::Bool

  function MPolyQuoRingElem(f::S, P::MPolyQuoRing{S}, simplified = false) where {S}
    @assert parent(f) === base_ring(P)
    return new{S}(f, P, simplified)
  end
end

@enable_all_show_via_expressify MPolyQuoRingElem

function AbstractAlgebra.expressify(a::MPolyQuoRingElem; context = nothing)
  return expressify(a.f, context = context)
end

function Base.deepcopy_internal(a::MPolyQuoRingElem, dict::IdDict)
  return MPolyQuoRingElem(Base.deepcopy_internal(a.f, dict), a.P, a.simplified)
end

########################################################################
# Representatives of elements in quotient rings and normal forms 
#
# Elements [a] ∈ P/I in quotients of polynomial rings 
# P = 𝕜[x₁,…,xₙ] by ideals I = ⟨f₁,…,fᵣ⟩ admit unique representatives 
# a ∈ P whenever a normal form algorithm exists for the polynomial 
# ring P. This is really a question about the ring of coefficients 𝕜. 
#
# In general, we can not expect a normal form algorithm to exist and, 
# in particular, that it is implemented in Singular. However, we wish 
# to use Singular as a default backend and this also drives the design 
# of the quotient rings to begin with. 
#
# To make sure that the data structure for quotient rings can 
# nevertheless also accommodate more exotic coefficient rings we
# provide the following functionality to decide the existence and use 
# of a Singular backend depending on the type.
#
# Since the concept of traits in Julia is not uniform, we briefly 
# describe how to set up your own type of coefficients with this
# framework. 
#
# Say you have a new type `MyType` for the `coefficient_ring` of a 
# polynomial ring `P` and you would like to make the `MPolyQuoRing`
# structure useful for quotients of the form `P/I`. Then you would 
# declare 
#
#   HasGroebnerAlgorithmTrait(::Type{MyType}) = HasSingularGroebnerAlgorithm()
#
# in case you are sure that the generic code to use Singular as a 
# backend can also digest polynomial rings whose coefficient ring 
# is of type `MyType`. 
#
# If you need to implement your own backend, you do the following.
# You declare
#
#   HasGroebnerAlgorithmTrait(::Type{MyType}) = HasMyCustomBackend()
#
# where `HasMyCustomBackend` is a name of your choice and then implement
#
# function _simplify(::HasMyCustomBackend, f::MPolyQuoRingElem)
#   # Do whatever has to be done to achieve a unique representative
#   ...
# end
#
# function _hash(::HasMyCustomBackend, f::MPolyQuoRingElem, u::UInt)
#   # Implement a hash which is unique for the class, NOT the 
#   # representative!
#   ...
# end
#
# function _is_equal(::HasMyCustomBackend, f::MPolyRingElem, g::MPolyQuoRingElem)
#   # Implement an equality check.
#   ...
# end
#
########################################################################

# The trait of the coefficient ring and its elements to decide 
# which backend to use for Groebner basis and normal form computations.
abstract type HasGroebnerAlgorithmTrait end
# Having a normal form algorithm available is strictly stronger than 
# having a Grobner basis algorithm. Thus we derive one type from the
# other here. 
abstract type HasNormalFormTrait <: HasGroebnerAlgorithmTrait end

# A normal form algorithm (for possibly non-global orderings) requires 
# strictly more, so we make it a subcase.
struct HasSingularNormalForm <: HasNormalFormTrait end
struct HasNoNormalForm <: HasNormalFormTrait end

struct HasSingularGroebnerAlgorithm <: HasGroebnerAlgorithmTrait end
struct HasRingFlattening <: HasGroebnerAlgorithmTrait end
struct HasNoGroebnerAlgorithm <: HasGroebnerAlgorithmTrait end

# By default we do not expect the Singular backend to be able to compute normal forms
HasNormalFormTrait(::Type{T}) where {T} = HasNoNormalForm()

# for convenience, allow passing in rings, ring elements or ring element types
HasNormalFormTrait(a::Ring) = HasNormalFormTrait(typeof(a))
HasNormalFormTrait(a::RingElem) = HasNormalFormTrait(parent_type(a))
HasNormalFormTrait(::Type{T}) where {T <: RingElem} = HasNormalFormTrait(parent_type(T))

# Same for the Groebner bases. But if a normal form algorithm exists, then 
# Groebner bases can be computed, too.
HasGroebnerAlgorithmTrait(a::Ring) = HasGroebnerAlgorithmTrait(typeof(a))
HasGroebnerAlgorithmTrait(a::RingElem) = HasGroebnerAlgorithmTrait(parent_type(a))
HasGroebnerAlgorithmTrait(::Type{T}) where {T <: RingElem} = HasGroebnerAlgorithmTrait(parent_type(T))

function HasGroebnerAlgorithmTrait(::Type{T}) where {T}
  HasNormalFormTrait(T) isa HasSingularNormalForm && return HasSingularGroebnerAlgorithm()
  return HasNoGroebnerAlgorithm()
end

# For polynomial rings over fields we expect Singular to be able to compute normal forms
HasNormalFormTrait(::Type{<:Field}) = HasSingularNormalForm()

# For polynomial rings over the integers we can compute Gröbner bases but
# in general can not compute normal forms.
HasGroebnerAlgorithmTrait(::Type{ZZRing}) = HasSingularGroebnerAlgorithm()
HasGroebnerAlgorithmTrait(::Type{zzModRing}) = HasSingularGroebnerAlgorithm()

# This list can (and should) be extended by eventual new types which 
# are supposed to make use of the Singular backend; see above.

# In particular, it is decided based on this trait whether a reasonable 
# hash function for elements in the quotient ring exists.

# In some cases it is useful to know that we can do a RingFlattening 
# and carry out certain operations in the de-nested ring
# These declarations happen in the file src/Rings/MPolyMap/flattenings.jl
# and need to be postponed due to limitations in inclusion orders.

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
  dim::Union{Int, Nothing, NegInf}
  gb::IdealGens{T}
  qRing::MPolyQuoRing

  function MPolyQuoIdeal(Ox::MPolyQuoRing{T}, si::Singular.sideal) where T <: MPolyRingElem
   @req singular_quotient_ring(Ox) == base_ring(si) "base rings must match"
   r = new{T}()
   r.gens  = IdealGens(Ox, si)
   r.qRing = Ox
   r.dim   = nothing
   R = base_ring(Ox)
   r.gens.O = [R(g) for g = gens(r.gens.S)]
   B = r.gens
   if length(B) >= 1 && is_graded(R)
     @req all(is_homogeneous, B.gens.O) "The generators of the ideal must be homogeneous"
   end
   return r
  end

  function MPolyQuoIdeal(Ox::MPolyQuoRing{T}, I::MPolyIdeal{T}) where T <: MPolyRingElem
    @req base_ring(Ox) === base_ring(I) "base rings must match"
    r = new{T}()
    r.gens = IdealGens(Ox, gens(I))
    r.qRing = Ox
    r.dim = nothing
    return r
  end
  
  function MPolyQuoIdeal(Ox::MPolyQuoRing{T}, V::Vector{T}) where T <: MPolyRingElem
    R = base_ring(Ox)
    @assert all(x->parent(x) == R, V)
    return MPolyQuoIdeal(Ox, ideal(R, V))
  end
end

function AbstractAlgebra.expressify(a::MPolyQuoIdeal; context = nothing)
  return Expr(:call, :ideal, [expressify(g, context = context) for g in gens(a)]...)
end

@doc raw"""
    base_ring(a::MPolyQuoIdeal)

Return the ambient ring of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> Q, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> a = ideal(Q, [x, y])
Ideal generated by
  x
  y

julia> base_ring(a)
Quotient
  of multivariate polynomial ring in 3 variables x, y, z
    over rational field
  by ideal (-x^2 + y, -x^3 + z)
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

function _groebner_basis(a::MPolyQuoIdeal)
  # Make sure that a has a Groebner basis?
  if !isdefined(a, :gb)
    a.gb = IdealGens(base_ring(a), Singular.std(singular_generators(a.gens)))
    a.gb.gens.S.isGB = a.gb.isGB = true
  end
end

function singular_groebner_generators(a::MPolyQuoIdeal)
  _groebner_basis(a)

  return a.gb.S
end

@doc raw"""
    gens(a::MPolyQuoIdeal)

Return the generators of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]))
(Quotient of multivariate polynomial ring by ideal (-x^2 + y, -x^3 + z), Map: R -> A)

julia> a = ideal(A, [x-y])
Ideal generated by
  x - y

julia> gens(a)
1-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x - y
```
"""
function gens(a::MPolyQuoIdeal)
  oscar_assure(a)
  return map(a.gens.Ox, a.gens.O)
end

gen(a::MPolyQuoIdeal, i::Int) = a.gens.Ox(a.gens.O[i])

@doc raw"""
    number_of_generators(a::MPolyQuoIdeal)

Return the number of generators of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]))
(Quotient of multivariate polynomial ring by ideal (-x^2 + y, -x^3 + z), Map: R -> A)

julia> a = ideal(A, [x-y])
Ideal generated by
  x - y

julia> number_of_generators(a)
1
```
"""
function number_of_generators(a::MPolyQuoIdeal)
  oscar_assure(a)
  return length(a.gens.O)
end


# powers, addition and multiplication do not require the singular quotient ring

@doc raw"""
    :^(a::MPolyQuoIdeal, m::Int)

Return the `m`-th power of `a`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, [x^2-y, y^2-x+y]);

julia> a = ideal(A, [x+y])
Ideal generated by
  x + y

julia> a^2
Ideal generated by
  x^2 + 2*x*y + y^2
```
"""
function Base.:^(a::MPolyQuoIdeal, m::Int)
  return MPolyQuoIdeal(base_ring(a), singular_generators(a.gens)^m)
end

@doc raw"""
    :+(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return the sum of `a` and `b`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, [x^2-y, y^2-x+y]);

julia> a = ideal(A, [x+y])
Ideal generated by
  x + y

julia> b = ideal(A, [x^2+y^2, x+y])
Ideal generated by
  x^2 + y^2
  x + y

julia> a+b
Ideal generated by
  x + y
  x^2 + y^2
```
"""
function Base.:+(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  @req base_ring(a) == base_ring(b) "base rings must match"
  return MPolyQuoIdeal(base_ring(a), singular_generators(a.gens) + singular_generators(b.gens))
end

@doc raw"""
    :*(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return the product of `a` and `b`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, [x^2-y, y^2-x+y]);

julia> a = ideal(A, [x+y])
Ideal generated by
  x + y

julia> b = ideal(A, [x^2+y^2, x+y])
Ideal generated by
  x^2 + y^2
  x + y

julia> a*b
Ideal generated by
  x^3 + x^2*y + x*y^2 + y^3
  x^2 + 2*x*y + y^2
```
"""
function Base.:*(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  @req base_ring(a) == base_ring(b) "base rings must match"
  return MPolyQuoIdeal(base_ring(a), singular_generators(a.gens) * singular_generators(b.gens))
end

@doc raw"""
    intersect(a::MPolyQuoIdeal{T}, bs::MPolyQuoIdeal{T}...) where T
    intersect(V::Vector{MPolyQuoIdeal{T}}) where T

Return the intersection of two or more ideals.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, ideal(R, [x^2-y^3, x-y]));

julia> a = ideal(A, [y^2])
Ideal generated by
  y^2

julia> b = ideal(A, [x])
Ideal generated by
  x

julia> intersect(a,b)
Ideal generated by
  x*y

julia> intersect([a,b])
Ideal generated by
  x*y
```
"""
function intersect(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}...) where T
  for g in b
    @req base_ring(g) == base_ring(a) "base rings must match"
  end
  as = Singular.intersection(singular_generators(a.gens), [singular_generators(g.gens) for g in b]...)
  return MPolyQuoIdeal(base_ring(a), as)
end

function intersect(V::Vector{MPolyQuoIdeal{T}}) where T
  @assert length(V) != 0
  length(V) == 1 && return V[1]

  return intersect(V[1], V[2:end]...)
end

#######################################################

@doc raw"""
    quotient(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return the ideal quotient of `a` by `b`. Alternatively, use `a:b`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, ideal(R, [x^2-y^3, x-y]));

julia> a = ideal(A, [y^2])
Ideal generated by
  y^2

julia> b = ideal(A, [x])
Ideal generated by
  x

julia> a:b
Ideal generated by
  y
```
"""
function quotient(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  @req base_ring(a) == base_ring(b) "base rings must match"
  return MPolyQuoIdeal(base_ring(a), Singular.quotient(singular_generators(a.gens), singular_generators(b.gens)))
end
(::Colon)(a::MPolyQuoIdeal, b::MPolyQuoIdeal) = quotient(a, b)

# TODO: replace by a more efficient method!
@attr Bool function is_prime(I::MPolyQuoIdeal)
  return is_prime(saturated_ideal(I))
end

@attr typeof(I) function radical(                                                   
    I::MPolyQuoIdeal; 
    preprocessing::Bool=true,
    eliminate_variables::Bool=true, 
    factor_generators::Bool=true
  )
  R = base_ring(I)
  J = saturated_ideal(I)
  return ideal(R, [g for g in R.(gens(radical(J; preprocessing, eliminate_variables, factor_generators))) if !iszero(g)])
end

# The following is to streamline the programmer's
# interface for the use of the four standard rings
# for the schemes `MPolyRing`, `MPolyQuoRing`, `MPolyLocRing`,
# and `MPolyQuoLocRing` together with their ideals.
# We return the preimage of the given ideal under the
# canonical map from the underlying free polynomial ring.
@attr Any function saturated_ideal(I::MPolyQuoIdeal)
  R = base_ring(base_ring(I))
  J = ideal(R, lift.(gens(I))) + modulus(base_ring(I))
  return J
end

# TODO: Replace by a more efficient method!
@attr Bool function is_radical(I::MPolyQuoIdeal)
  return is_radical(saturated_ideal(I))
end

@doc raw"""
    is_zero(a::MPolyQuoIdeal)

Return `true` if `a` is the zero ideal, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, [x^2-y, y^2-x+y]);

julia> a = ideal(A, [x^2+y^2, x+y])
Ideal generated by
  x^2 + y^2
  x + y

julia> is_zero(a)
false

julia> b = ideal(A, [x^2-y])
Ideal generated by
  x^2 - y

julia> is_zero(b)
true
```
"""
@attr Bool function is_zero(a::MPolyQuoIdeal)
  R = base_ring(a)
  return Singular.iszero(Singular.reduce(singular_generators(a.gens), singular_quotient_groebner_basis(R)))
end

@doc raw"""
    ideal(A::MPolyQuoRing{T}, V::Vector{T}) where T <: MPolyRingElem

Given a (graded) quotient ring `A=R/I` and a vector `V` of (homogeneous) polynomials in `R`,
create the ideal of `A` which is generated by the images of the entries of `V`.

    ideal(A::MPolyQuoRing{T}, V::Vector{MPolyQuoRingElem{T}}) where T <: MPolyRingElem

Given a (graded) quotient ring `A` and a vector `V` of (homogeneous) elements of `A`,
create the ideal of `A` which is generated by the entries of `V`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, ideal(R, [x^2-y^3, x-y]));

julia> I = ideal(A, [x^2-y])
Ideal generated by
  x^2 - y
```
```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> B, _ = quo(S, ideal(S, [x^2*z-y^3, x-y]));

julia> J = ideal(B, [x^2-y^2])
Ideal generated by
  x^2 - y^2
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

ideal(R::MPolyQuoRing, V::Vector) = ideal(R, elem_type(R)[R(x) for x in V])


function ideal(A::MPolyQuoRing{T}, x::T) where T <: MPolyRingElem
  return ideal(A,[x])
end

function ideal(A::MPolyQuoRing{T}, x::MPolyQuoRingElem{T}) where T <: MPolyRingElem
  return ideal(A,[x])
end
##################################################################

parent_type(::Type{MPolyQuoRingElem{S}}) where S = MPolyQuoRing{S}
elem_type(::Type{MPolyQuoRing{S}})  where S= MPolyQuoRingElem{S}

canonical_unit(a::MPolyQuoRingElem) = one(parent(a))

parent(a::MPolyQuoRingElem) = a.P

zero(R::MPolyQuoRing) = MPolyQuoRingElem(zero(base_ring(R)), R, true)

function is_zero(a::MPolyQuoRingElem)
  return iszero(simplify(a).f)
end


function check_parent(a::MPolyQuoRingElem, b::MPolyQuoRingElem)
  @req parent(a) == parent(b) "parents must match"
  return true
end

+(a::MPolyQuoRingElem{S}, b::MPolyQuoRingElem{S}) where {S} = check_parent(a, b) && MPolyQuoRingElem(a.f+b.f, a.P)

-(a::MPolyQuoRingElem, b::MPolyQuoRingElem) = check_parent(a, b) && MPolyQuoRingElem(a.f-b.f, a.P)

-(a::MPolyQuoRingElem) = MPolyQuoRingElem(-a.f, a.P)

*(a::MPolyQuoRingElem{S}, b::MPolyQuoRingElem{S}) where {S} = check_parent(a, b) && simplify(MPolyQuoRingElem(a.f*b.f, a.P))

function Base.:(^)(a::MPolyQuoRingElem, b::Base.Integer)
  if b >= 0
    simplify(MPolyQuoRingElem(Base.power_by_squaring(a.f, b), a.P))
  else
    return inv(a)^(-b)
  end
end

*(a::MPolyQuoRingElem, b::QQFieldElem) = simplify(MPolyQuoRingElem(a.f * b, a.P))

*(a::MPolyQuoRingElem, b::ZZRingElem) = simplify(MPolyQuoRingElem(a.f * b, a.P))

*(a::QQFieldElem, b::MPolyQuoRingElem) = simplify(MPolyQuoRingElem(a * b.f, b.P))

*(a::ZZRingElem, b::MPolyQuoRingElem) = simplify(MPolyQuoRingElem(a * b.f, b.P))

#*(a::MPolyQuoRingElem, b::MPolyQuoRingElem) = check_parent(a, b) && MPolyQuoRingElem(a.f*b.f, a.P)
#
#^(a::MPolyQuoRingElem, b::Base.Integer) = MPolyQuoRingElem(Base.power_by_squaring(a.f, b), a.P)

function Oscar.mul!(a::MPolyQuoRingElem, b::MPolyQuoRingElem, c::MPolyQuoRingElem)
  a.f = b.f * c.f
  a.simplified = false
  return a
end

function Oscar.add!(a::MPolyQuoRingElem, b::MPolyQuoRingElem, c::MPolyQuoRingElem)
  a.f = b.f + c.f
  a.simplified = false
  return a
end

@doc raw"""
    simplify(a::MPolyQuoIdeal)

If `a` is an ideal of the affine algebra `A = R/I`, say, replace the internal polynomial representative 
of each generator of `a` by its normal form mod `I` with respect to `ordering(A)`.


# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> a = ideal(A, [x^3*y^4-x+y, x*y+y^2*x])
Ideal generated by
  x^3*y^4 - x + y
  x*y^2 + x*y

julia> gens(a)
2-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x^3*y^4 - x + y
 x*y^2 + x*y

julia> simplify(a)
Ideal generated by
  x^2*y^3 - x + y
  x*y^2 + x*y

julia> gens(a)
2-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x^2*y^3 - x + y
 x*y^2 + x*y
```
"""
function simplify(a::MPolyQuoIdeal)
  Q = base_ring(a)
  R = base_ring(Q)
  red  = reduce(singular_generators(a.gens), singular_quotient_groebner_basis(Q))
  SQ   = singular_poly_ring(Q)
  si   = Singular.Ideal(SQ, unique!(gens(red)))
  a.gens.S = si
  a.gens.O = [R(g) for g = gens(a.gens.S)]
  return a
end

#######################################################

@doc raw"""
    ideal_membership(f::MPolyQuoRingElem{T}, a::MPolyQuoIdeal{T}) where T

Return `true` if `f` is contained in `a`, `false` otherwise. Alternatively, use `f in a`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> a = ideal(A, [x^3*y^4-x+y, x*y+y^2*x])
Ideal generated by
  x^3*y^4 - x + y
  x*y^2 + x*y

julia> f = A(x^2*y^3-x+y)
x^2*y^3 - x + y

julia> f in a
true
```
"""
function ideal_membership(a::MPolyQuoRingElem{T}, b::MPolyQuoIdeal{T}) where T
  parent(a) == base_ring(b) || error("base rings must match")
  SR = singular_poly_ring(base_ring(b))
  return Singular.iszero(Singular.reduce(SR(simplify(a)), singular_groebner_generators(b)))
end

Base.:in(a::MPolyQuoRingElem, b::MPolyQuoIdeal) = ideal_membership(a, b)


@doc raw"""
    is_subset(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return `true` if `a` is contained in `b`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> a = ideal(A, [x^3*y^4-x+y, x*y+y^2*x])
Ideal generated by
  x^3*y^4 - x + y
  x*y^2 + x*y

julia> b = ideal(A, [x^3*y^3-x+y, x^2*y+y^2*x])
Ideal generated by
  x^3*y^3 - x + y
  x^2*y + x*y^2

julia> is_subset(a,b)
false

julia> is_subset(b,a)
true
```
"""
function is_subset(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  @req base_ring(a) == base_ring(b) "base rings must match"
  simplify(a)
  return Singular.iszero(Singular.reduce(singular_generators(a.gens), singular_groebner_generators(b)))
end

@doc raw"""
    ==(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T

Return `true` if `a` is equal to `b`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, _ = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> a = ideal(A, [x^3*y^4-x+y, x*y+y^2*x])
Ideal generated by
  x^3*y^4 - x + y
  x*y^2 + x*y

julia> b = ideal(A, [x^3*y^3-x+y, x^2*y+y^2*x])
Ideal generated by
  x^3*y^3 - x + y
  x^2*y + x*y^2

julia> a == b
false
```
"""
function Base.:(==)(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
  return issubset(a, b) && issubset(b, a)
end

@doc raw"""
    simplify(f::MPolyQuoRingElem)

If `f` is an element of the affine algebra `A = R/I`, say, replace the internal polynomial representative of `f` 
by its normal form mod `I` with respect to `ordering(A)`.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, [:x]);

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
  Q = parent(f)::MPolyQuoRing
  P = base_ring(Q)::MPolyRing
  kk = coefficient_ring(P)::Ring
  return _simplify(HasGroebnerAlgorithmTrait(kk), f)
end

function _simplify(::HasSingularGroebnerAlgorithm, f::MPolyQuoRingElem)
  f.simplified && return f
  R  = parent(f)
  OR = oscar_origin_ring(R)
  SR = singular_origin_ring(R)
  G  = singular_origin_groebner_basis(R)
  g  = f.f
  f.f = OR(reduce(SR(g), G))
  f.simplified = true
  return f::elem_type(R)
end

# The above method for `simplify` assume that there is a singular backend which 
# can be used. However, we are using (graded) quotient rings also with coefficient 
# rings R which can not be translated to Singular; for instance when R is again 
# a polynomial ring, or a quotient/localization thereof, or even an `AffineSchemeOpenSubschemeRing`. 
# Still in many of those cases, we can use `RingFlattening` to bind a computational 
# backend. In particular, this allows us to do ideal_membership tests; see 
# the file `flattenings.jl` for details. 
#
# The generic method below is a compromise in the sense that `simplify` does not reduce 
# a given element to a unique representative as would be the case in a groebner basis reduction, 
# but it nevertheless reduces the element to zero in case its representative is 
# contained in the modulus. This allows for both, the use of `RingFlattening`s and 
# the potential speedup of `iszero` tests. 
function _simplify(::HasRingFlattening, f::MPolyQuoRingElem)
  f.simplified && return f
  if f.f in modulus(parent(f))
    f.f = zero(f.f)
  end
  f.simplified = true
  return f::elem_type(parent(f))
end

# Having the simplify method do nothing still allows for arithmetic to operate 
# (simplify is automatically called for multiplications), but further functionality 
# like equality and hashing will throw an error.
function _simplify(::HasNoGroebnerAlgorithm, f::MPolyQuoRingElem)
  return f
  # error("no groebner backend available for simplification; you can specify a groebner backend by implementing the `HasGroebnerAlgorithmTrait(::Type{T}) = MyBackend()` for `T = $(typeof(coefficient_ring(base_ring(parent(f)))))`")
end

@doc raw"""
    ==(f::MPolyQuoRingElem{T}, g::MPolyQuoRingElem{T}) where T

Return `true` if `f` is equal to `g`, `false` otherwise.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, [:x]);

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
  Q = parent(f)::MPolyQuoRing
  P = base_ring(Q)::MPolyRing
  kk = coefficient_ring(P)::Ring
  return _is_equal(HasGroebnerAlgorithmTrait(kk), f, g)
end

function _is_equal(::HasSingularGroebnerAlgorithm, f::MPolyQuoRingElem{T}, g::MPolyQuoRingElem{T}) where T
  f.f == g.f && return true
  simplify(f)
  simplify(g)
  return f.f == g.f
end

# By default we refer to the generic ideal membership routine which 
# might be implemented by other means, for instance via a `RingFlattening`.
function _is_equal(::HasRingFlattening, f::MPolyQuoRingElem{T}, g::MPolyQuoRingElem{T}) where T
  return f.f - g.f in modulus(parent(f))
end

function _is_equal(::HasNoGroebnerAlgorithm, f::MPolyQuoRingElem{T}, g::MPolyQuoRingElem{T}) where T
  error("no groebner backend available for equality check; you can specify a groebner backend by implementing the `HasGroebnerAlgorithmTrait(::Type{T}) = MyBackend()` for `T = $(typeof(coefficient_ring(base_ring(parent(f)))))`")
end

@doc raw"""
    quo(R::MPolyRing, I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(R)) -> MPolyQuoRing, Map

Create the quotient ring `R/I` and return the new ring as well as the projection map `R` $\to$ `R/I`.

    quo(R::MPolyRing, V::Vector{MPolyRingElem}; ordering::MonomialOrdering = default_ordering(R)) -> MPolyQuoRing, Map

As above, where `I` is the ideal of `R` generated by the polynomials in `V`.

!!! note
    Once `R/I` is created,  all computations within `R/I` relying on division with remainder and/or Gröbner bases are done with respect to `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> A, p = quo(R, ideal(R, [x^2-y^3, x-y]));

julia> A
Quotient
  of multivariate polynomial ring in 2 variables x, y
    over rational field
  by ideal (x^2 - y^3, x - y)

julia> typeof(A)
MPolyQuoRing{QQMPolyRingElem}

julia> typeof(x)
QQMPolyRingElem

julia> p
Map defined by a julia-function with inverse
  from multivariate polynomial ring in 2 variables over QQ
  to quotient of multivariate polynomial ring by ideal (x^2 - y^3, x - y)

julia> p(x)
x

julia> typeof(p(x))
MPolyQuoRingElem{QQMPolyRingElem}
```
```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> B, _ = quo(S, ideal(S, [x^2*z-y^3, x-y]))
(Quotient of multivariate polynomial ring by ideal (x^2*z - y^3, x - y), Map: S -> B)

julia> typeof(B)
MPolyQuoRing{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}
```
"""
function quo(R::MPolyRing, I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(R))
  q = MPolyQuoRing(R, I, ordering)
  function im(a::MPolyRingElem)
    parent(a) !== R && error("Element not in the domain of the map")
    return MPolyQuoRingElem(a, q)
  end
  function pr(a::MPolyQuoRingElem)
    return a.f
  end
  return q, MapFromFunc(R, q, im, pr)
end

function quo(R::MPolyRing, I::Vector{<:MPolyRingElem}; ordering::MonomialOrdering = default_ordering(R))
  return quo(R, ideal(I); ordering)
end

function quo(R::MPolyRing, f::MPolyRingElem...; ordering::MonomialOrdering = default_ordering(R))
  return quo(R, ideal(collect(f)); ordering)
end

lift(a::MPolyQuoRingElem) = a.f


(Q::MPolyQuoRing)() = MPolyQuoRingElem(base_ring(Q)(), Q)

function (Q::MPolyQuoRing)(a::MPolyQuoRingElem; check::Bool=true)
  if parent(a) === Q
    return a
  else
    return Q(base_ring(Q)(a))
  end
end

function(Q::MPolyRing{T})(a::MPolyQuoRingElem{<:MPolyRingElem{T}}; check::Bool=false) where {T}
  @req base_ring(parent(a)) === Q "parent mismatch"
  return lift(a)
end


function (Q::MPolyQuoRing{S})(a::S; check::Bool=false) where {S <: MPolyRingElem}
  @req base_ring(Q) === parent(a) "parent mismatch"
  return MPolyQuoRingElem(a, Q)
end

function (Q::MPolyQuoRing)(a::MPolyRingElem; check::Bool=false)
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

(Q::MPolyQuoRing)(a; check::Bool=false) = MPolyQuoRingElem(base_ring(Q)(a), Q)

one(Q::MPolyQuoRing) = Q(1)

@doc raw"""
    is_invertible_with_inverse(f::MPolyQuoRingElem)

If `f` is invertible with inverse `g`, say, return `(true, g)`. Otherwise, return `(false, f)`.

# Examples

```jldoctest
julia> R, c = polynomial_ring(QQ, :c => (1:3));

julia> R, c = grade(R, [1, 2, 3]);

julia> I = ideal(R, [  -c[1]^3 + 2*c[1]*c[2] - c[3], c[1]^4 - 3*c[1]^2*c[2] + 2*c[1]*c[3] + c[2]^2,-c[1]^5 + 4*c[1]^3*c[2] - 3*c[1]^2*c[3] - 3*c[1]*c[2]^2 + 2*c[2]*c[3]]);

julia> A, _ = quo(R, I);

julia> f = A(c[1]^2 - c[1] - c[2] + 1)
c[1]^2 - c[1] - c[2] + 1

julia> tt, g = is_invertible_with_inverse(f)
(true, c[1] + c[2] + c[3] + 1)

julia> f*g
1

```
"""
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

  if Q isa MPolyQuoRing{<:MPolyDecRingElem}
     J = [x.f for x in J]
  end

  j, T = standard_basis_with_transformation_matrix(ideal(J))
  if is_constant(j[1]) && is_unit(first(coefficients(j[1])))
    @assert ncols(T) == 1
    return true, inv(first(coefficients(j[1])))*Q(T[end, 1])
  end
  return false, a
end

is_unit(a::MPolyQuoRingElem) = is_invertible_with_inverse(a)[1]

@doc raw"""
    inv(f::MPolyQuoRingElem)

If `f` is invertible, return its inverse. Otherwise, throw an error.

# Examples

```jldoctest
julia> R, c = polynomial_ring(QQ, :c => (1:3));

julia> I = ideal(R, [  -c[1]^3 + 2*c[1]*c[2] - c[3], c[1]^4 - 3*c[1]^2*c[2] + 2*c[1]*c[3] + c[2]^2,-c[1]^5 + 4*c[1]^3*c[2] - 3*c[1]^2*c[3] - 3*c[1]*c[2]^2 + 2*c[2]*c[3]]);

julia> A, _ = quo(R, I);

julia> f = A(c[1]^2 - c[1] - c[2] + 1)
c[1]^2 - c[1] - c[2] + 1

julia> inv(f)
c[1] + c[2] + c[3] + 1

```
"""
function inv(a::MPolyQuoRingElem)
  fl, b = is_invertible_with_inverse(a)
  fl || error("Element not invertible")
  return b
end

### The following method is only implemented when the coefficient ring is a field.
# The code should be valid generically, but the Singular backend needed for the
# ideal quotient is probably buggy for non-fields.
function is_zero_divisor(f::MPolyQuoRingElem{<:MPolyRingElem{<:FieldElem}})
  iszero(f) && return true
  b = simplify(f)
  # The next block is basically useless when the coefficient ring is
  # a field, because it is merely another `is_zero`-check. However,
  # once more functionality is working, it will actually do stuff and
  # the above signature can be widened.
  if is_constant(lift(b))
    c = first(AbstractAlgebra.coefficients(lift(b)))
    return is_zero_divisor(c)
  end
  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
end

"""
Converts a sparse-Singular vector of polynomials to an Oscar sparse row.
"""
function sparse_row(R::MPolyRing, M::Singular.svector{<:Singular.spoly})
  v = Dict{Int, MPolyBuildCtx}()
  for (i, e, c) = M
    vi = get!(v, i) do
      MPolyBuildCtx(R)
    end
    push_term!(vi, base_ring(R)(c), e)
  end
  pos_value_vector::Vector{Tuple{Int, elem_type(R)}} = [(k,finish(v)) for (k,v) = v]
  return sparse_row(R, pos_value_vector)
end

"""
Converts a sparse-Singular vector of polynomials to an Oscar sparse row.
Collect only the column indices in `U`.
"""
function sparse_row(R::MPolyRing, M::Singular.svector{<:Singular.spoly}, U::AbstractUnitRange)
  v = Dict{Int, MPolyBuildCtx}()
  for (i, e, c) = M
    (i in U) || continue
    vi = get!(v, i) do
      MPolyBuildCtx(R)
    end
    push_term!(vi, base_ring(R)(c), e)
  end
  pos_value_vector::Vector{Tuple{Int, elem_type(R)}} = [(k,finish(v)) for (k,v) = v]
  return sparse_row(R, pos_value_vector)
end

"""
Converts the sparse-Singular matrix (`Module`) row by row to an Oscar sparse-matrix.
Only the row indices (generators) in `V` and the column indices in `U` are converted.
"""
function sparse_matrix(R::MPolyRing, M::Singular.Module, V::AbstractUnitRange, U::AbstractUnitRange)
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

  J = vcat(J, [b.f])
  BJ = IdealGens(J, keep_ordering = false)

  s, rest = Singular.lift(singular_generators(BJ), singular_generators(BS))
  if !iszero(rest)
    return false, a
  end
  return true, Q(sparse_matrix(base_ring(Q), s, 1:1, length(J):length(J))[1, length(J)])
end

function divexact(a::MPolyQuoRingElem, b::MPolyQuoRingElem; check::Bool=true)
  check_parent(a, b)
  b, q = divides(a, b)
  !b && error("Division is not exact in divexact")
  return q
end

### 
# The following two functions below provide a hotfix to make sure that the preferred 
# ordering provided to the constructor of the quotient ring is actually used for the 
# groebner basis computations.
function ordering(A::MPolyQuoRing)
  return A.ordering
end

function _divides_hack(a::MPolyQuoRingElem, b::MPolyQuoRingElem)
  check_parent(a, b)
  iszero(a) && iszero(b) && return (true, zero(parent(a)))
  iszero(b) && error("cannot divide by zero")

  A = parent(a)
  R = base_ring(A)
  mod = modulus(A)
  o = ordering(A)
  # Take a groebner basis for the preferred ordering, hoping that it will speed up things
  gb_mod = gens(groebner_basis(mod, ordering=o))
  I = ideal(R, push!(gb_mod, lift(b)))
  # Make sure that the singular side of this ideal is filled with the correct ordering
  # Get our hands on the actual singular ideal
  Ising = singular_generators(I.gens, o)
  # ...and the ring
  Rsing = base_ring(Ising)
  a_ideal = Singular.Ideal(Rsing, [Rsing(lift(a))])
  u_sing, rem = Singular.lift(Ising, a_ideal)
  !iszero(rem) && return false, a
  return true, A(sparse_matrix(base_ring(A), u_sing, 1:1, ngens(I):ngens(I))[1, ngens(I)])
end

#TODO: find a more descriptive, meaningful name
function _kbase(Q::MPolyQuoRing)
  G = singular_origin_groebner_basis(Q)
  s = Singular.kbase(G)
  if iszero(s)
    error("ideal is not zero-dimensional")
  end
  return [base_ring(Q)(x) for x = gens(s)]
end

function vector_space(K::AbstractAlgebra.Field, Q::MPolyQuoRing)
  R = base_ring(Q)
  @assert K == base_ring(R)
  l = _kbase(Q)
  V = free_module(K, length(l))
  function im(a::Generic.FreeModuleElem)
    @assert parent(a) === V
    b = R(0)
    for k=1:length(l)
      c = a[k]
      if !iszero(c)
        b += c*l[k]
      end
    end
    return Q(b)
  end

  # The inverse function. We use the fact that for a chosen monomial ordering 
  # the monomials which are not in the leading ideal, form a basis for the 
  # quotient; see Greuel/Pfister "A singular introduction to Commutative Algebra".
  function prim(a::MPolyQuoRingElem)
    @assert parent(a) === Q
    b = lift(a)::MPolyRingElem
    o = default_ordering(R)
    # TODO: Make sure the ordering is the same as the one used for the _kbase above
    @assert is_global(o) "ordering must be global"
    b = normal_form(b, modulus(Q), ordering=o)
    result = zero(V)
    while !iszero(b)
      m = leading_monomial(b, ordering=o)
      c = leading_coefficient(b, ordering=o)
      j = findfirst(==(m), l)
      result = result + c * V[j]
      b = b - c * m
    end
    return result
  end
  return V, MapFromFunc(V, Q, im, prim)
end

# To fix printing of fraction fields of MPolyQuoRing
function AbstractAlgebra.expressify(a::AbstractAlgebra.Generic.FracFieldElem{T};
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

function _grading(R::MPolyQuoRing)
  if base_ring(R) isa MPolyDecRing
    return _grading(base_ring(R))
  else
    error("Underlying polynomial ring must be graded")
  end
end

@doc raw"""
    degree(f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given a homogeneous element `f` of a graded affine algebra, return the degree of `f`.

    degree(::Type{Vector{Int}}, f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given a homogeneous element `f` of a $\mathbb Z^m$-graded affine algebra, return the degree of `f`, converted to a vector of integer numbers.

    degree(::Type{Int}, f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given a homogeneous element `f` of a $\mathbb Z$-graded affine algebra, return the degree of `f`, converted to an integer number.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z] );

julia> A, p = quo(R, ideal(R, [y-x, z^3-x^3]))
(Quotient of multivariate polynomial ring by ideal (-x + y, -x^3 + z^3), Map: R -> A)

julia> f = p(y^2-x^2+z^4)
-x^2 + y^2 + z^4

julia> degree(f)
[4]

julia> typeof(degree(f))
FinGenAbGroupElem

julia> degree(Int, f)
4

julia> typeof(degree(Int, f))
Int64
```
"""
function degree(a::MPolyQuoRingElem{<:MPolyDecRingElem}; check::Bool=true)
  simplify(a)
  @req !iszero(a) "Element must be non-zero"
  return degree(a.f; check)
end

function _degree_fast(a::MPolyQuoRingElem{<:MPolyDecRingElem})
  simplify(a)
  @req !iszero(a) "Element must be non-zero"
  return _degree_fast(a.f)
end
function degree(::Type{Int}, a::MPolyQuoRingElem{<:MPolyDecRingElem}; check::Bool=true)
  @assert is_z_graded(base_ring(parent(a)))
  return Int(degree(a; check)[1])
end

function degree(::Type{Vector{Int}}, a::MPolyQuoRingElem{<:MPolyDecRingElem}; check::Bool=true)
  @assert is_zm_graded((base_ring(parent(a))))
  d = degree(a; check)
  return Int[d[i] for i=1:ngens(parent(d))]
end

is_filtered(q::MPolyQuoRing) = is_filtered(base_ring(q))
is_graded(q::MPolyQuoRing) = is_graded(base_ring(q))

@doc raw"""
    homogeneous_component(f::MPolyQuoRingElem{<:MPolyDecRingElem}, g::FinGenAbGroupElem)

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
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> A, p = quo(R, ideal(R, [y-x, z^3-x^3]));

julia> f = p(y^2-x^2+x*y*z+z^4)
-x^2 + x*y*z + y^2 + z^4

julia> homogeneous_component(f, 4)
z^4
```
"""
function homogeneous_component(a::MPolyQuoRingElem{<:MPolyDecRingElem}, d::FinGenAbGroupElem)
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

@doc raw"""
    homogeneous_components(f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given an element `f` of a graded affine algebra, return the homogeneous components of `f`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> A, p = quo(R, ideal(R, [y-x, z^3-x^3]));

julia> f = p(y^2-x^2+x*y*z+z^4)
-x^2 + x*y*z + y^2 + z^4

julia> homogeneous_components(f)
Dict{FinGenAbGroupElem, MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}} with 2 entries:
  [4] => z^4
  [3] => y^2*z
```
"""
function homogeneous_components(a::MPolyQuoRingElem{<:MPolyDecRingElem})
  simplify(a)
  h = homogeneous_components(a.f)
  return Dict{keytype(h), typeof(a)}(x => parent(a)(y) for (x, y) in h)
end

@doc raw"""
    is_homogeneous(f::MPolyQuoRingElem{<:MPolyDecRingElem})

Given an element `f` of a graded affine algebra, return `true` if `f` is homogeneous, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

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

@doc raw"""
    grading_group(A::MPolyQuoRing{<:MPolyDecRingElem})

If `A` is, say, `G`-graded, return `G`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> A, _ = quo(R, ideal(R, [x^2*z-y^3, x-y]));

julia> grading_group(A)
Z
```
"""
function grading_group(A::MPolyQuoRing{<:MPolyDecRingElem})
  return grading_group(base_ring(A))
end

function hash(w::MPolyQuoRingElem, u::UInt)
  Q = parent(w)::MPolyQuoRing
  P = base_ring(Q)::MPolyRing
  kk = coefficient_ring(P)::Ring
  return _hash(HasGroebnerAlgorithmTrait(kk), w, u)
end

function _hash(::HasSingularGroebnerAlgorithm, w::MPolyQuoRingElem, u::UInt)
  simplify(w)
  return hash(w.f, u)
end

function _hash(::HasGroebnerAlgorithmTrait, w::MPolyQuoRingElem, u::UInt)
  error("hash function not implemented due to lack of unique representatives; you can specify a groebner backend by implementing the `HasGroebnerAlgorithmTrait(::Type{T}) = MyBackend()` for `T = $(typeof(coefficient_ring(base_ring(parent(w)))))`")
end

################################################################
### homogeneous components quotient ring
################################################################

@doc raw"""
    homogeneous_component(A::MPolyQuoRing{<:MPolyDecRingElem}, g::FinGenAbGroupElem)

Given a graded quotient `A` of a multivariate polynomial ring over a field, 
where the grading group is free of type `FinGenAbGroup`, and given an element `g` of 
that group, return the homogeneous component of `A` of degree `g`. Additionally, return
the embedding of the component into `A`.

    homogeneous_component(A::MPolyQuoRing{<:MPolyDecRingElem}, g::Vector{<:IntegerUnion})

Given a $\mathbb  Z^m$-graded quotient `A` of a multivariate polynomial ring over a field, 
and given a vector `g` of $m$ integers, convert `g` into an element of the grading
group of `A`, and return the homogeneous component of `A` whose degree 
is that element. Additionally, return the embedding of the component into `A`.

    homogeneous_component(A::MPolyQuoRing{<:MPolyDecRingElem}, g::IntegerUnion)

Given a $\mathbb  Z$-graded quotient `A` of a multivariate polynomial ring over a field, 
and given an integer `g`, convert `g` into an element of the grading group of `A`, 
and return the homogeneous component of `A` whose degree is that element.
Additionally, return the embedding of the component into `A`.

!!! note
    If the component is not finite dimensional, an error message will be thrown.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z])
(Graded multivariate polynomial ring in 4 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w, x, y, z])

julia> L = homogeneous_component(R, 2);

julia> HC = gens(L[1]);

julia> EMB = L[2]
Map defined by a julia-function with inverse
  from R_[2] of dim 10
  to graded multivariate polynomial ring in 4 variables over QQ

julia> for i in 1:length(HC) println(EMB(HC[i])) end
z^2
y*z
y^2
x*z
x*y
x^2
w*z
w*y
w*x
w^2

julia> PTC = ideal(R, [-x*z + y^2, -w*z + x*y, -w*y + x^2]);

julia> A, _ = quo(R, PTC);

julia> L = homogeneous_component(A, 2);

julia> HC = gens(L[1]);

julia> EMB = L[2]
Map defined by a julia-function with inverse
  from quotient space over QQ with 7 generators and no relations
  to quotient of multivariate polynomial ring by ideal (-x*z + y^2, -w*z + x*y, -w*y + x^2)

julia> for i in 1:length(HC) println(EMB(HC[i])) end
z^2
y*z
x*z
w*z
w*y
w*x
w^2
```

```jldoctest
julia> G = abelian_group([0, 0])
Z^2

julia> W = [G[1], G[1], G[2], G[2], G[2]];

julia> S, x, y = graded_polynomial_ring(QQ, :x => 1:2, :y => 1:3; weights = W);

julia> L = homogeneous_component(S, [2,1]);

julia> HC = gens(L[1]);

julia> EMB = L[2]
Map defined by a julia-function with inverse
  from S_[2 1] of dim 9
  to graded multivariate polynomial ring in 5 variables over QQ

julia> for i in 1:length(HC) println(EMB(HC[i])) end
x[2]^2*y[3]
x[2]^2*y[2]
x[2]^2*y[1]
x[1]*x[2]*y[3]
x[1]*x[2]*y[2]
x[1]*x[2]*y[1]
x[1]^2*y[3]
x[1]^2*y[2]
x[1]^2*y[1]

julia> I = ideal(S, [x[1]*y[1]-x[2]*y[2]]);

julia> A, = quo(S, I);

julia> L = homogeneous_component(A, [2,1]);

julia> HC = gens(L[1]);

julia> EMB = L[2]
Map defined by a julia-function with inverse
  from quotient space over QQ with 7 generators and no relations
  to quotient of multivariate polynomial ring by ideal (x[1]*y[1] - x[2]*y[2])

julia> for i in 1:length(HC) println(EMB(HC[i])) end
x[2]^2*y[3]
x[2]^2*y[2]
x[2]^2*y[1]
x[1]*x[2]*y[3]
x[1]*x[2]*y[2]
x[1]^2*y[3]
x[1]^2*y[2]
```
"""
function homogeneous_component(W::MPolyQuoRing{<:MPolyDecRingElem}, d::FinGenAbGroupElem)
  #TODO: lazy: ie. no enumeration of points
  #      apparently it is possible to get the number of points faster than the points
  D = parent(d)
  @assert D == grading_group(W)
  R = base_ring(W)

  H, mH = homogeneous_component(R, d)
  I = modulus(W)
  M = gens(leading_ideal(I))
  cache = Dict{typeof(d), typeof(mH)}()
  q = Set{elem_type(H)}()
  for h = M
    g = degree(h)
    mI = get!(cache, g) do
      return homogeneous_component(R, d-g)[2]
    end
    for x = gens(domain(mI))
      push!(q, preimage(mH, h*mI(x)))
    end
  end

  s, ms = sub(H, collect(q))
  Q, mQ = quo(H, s)
#  set_attribute!(Q, :show => show_homo_comp, :data => (W, d))
  return Q, MapFromFunc(Q, W, x->W(mH((preimage(mQ, x)))), y->mQ(preimage(mH, y.f)))
end

function homogeneous_component(W::MPolyQuoRing{<:MPolyDecRingElem}, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(W)
  return homogeneous_component(W, grading_group(W)(g))
end

function homogeneous_component(W::MPolyQuoRing{<:MPolyDecRingElem}, g::IntegerUnion)
  @assert is_z_graded(W)
  return homogeneous_component(W, grading_group(W)([g]))
end

@doc raw"""
    dim(a::MPolyQuoIdeal)

Return the Krull dimension of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> a = ideal(A, [x-y])
Ideal generated by
  x - y

julia> dim(a)
0
```
"""
function dim(a::MPolyQuoIdeal)
  if a.dim === nothing 
    a.dim = Singular.dimension(singular_groebner_generators(a))
    a.dim == -1 && (a.dim = -inf)
  end
  return a.dim::Union{Int, NegInf}
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
@doc raw"""
    minimal_generating_set(I::MPolyQuoIdeal{<:MPolyDecRingElem})

Given a homogeneous ideal `I` of a graded affine algebra over a field,
return an array containing a minimal set of generators of `I`. If `I`
is the zero ideal, an empty list is returned.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> A, p = quo(R, ideal(R, [x-y]));

julia> V = [x, z^2, x^3+y^3, y^4, y*z^5];

julia> a = ideal(A, V);

julia> minimal_generating_set(a)
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 y
 z^2

julia> a = ideal(A, [x-y])
Ideal generated by
  x - y

julia> minimal_generating_set(a)
MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}[]
```
"""
function minimal_generating_set(I::MPolyQuoIdeal{<:MPolyDecRingElem})
  # This only works / makes sense for homogeneous ideals. So far ideals in an
  # MPolyDecRing are forced to be homogeneous though.
  Q = base_ring(I)

  @assert is_graded(Q)

  @req coefficient_ring(Q) isa AbstractAlgebra.Field "The coefficient ring must be a field"

  if isdefined(I, :gb)
    _, sing_min = Singular.mstd(singular_generators(I.gb.gens))
    return filter(!iszero, (Q).(gens(sing_min)))
  else
    sing_gb, sing_min = Singular.mstd(singular_generators(I.gens))
    I.gb = IdealGens(I.gens.Ox, sing_gb, true)
    I.gb.gens.S.isGB = I.gb.isGB = true
    return filter(!iszero, (Q).(gens(sing_min)))
  end
end

@doc raw"""
    small_generating_set(I::MPolyQuoIdeal)

Given a ideal `I` of an affine algebra over a field, return an array
containing set of generators of `I`, which is usually smaller than the
original one. If `I` is the zero ideal an empty list is returned.

!!! note
   Minimal generating sets exist only in the local and the homogeneous case. Beyond these cases, the best one can hope for is some small set of generators.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A, p = quo(R, ideal(R, [x-y]));

julia> V = [x, z^2, x^3+y^3, y^4, y*z^5];

julia> a = ideal(A, V);

julia> small_generating_set(a)
2-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 y
 z^2

```
"""
@attr Vector{elem_type(base_ring(I))} function small_generating_set(
    I::MPolyQuoIdeal; 
    algorithm::Symbol=:simple
  )
  # For non-homogeneous ideals, we do not have a notion of minimal generating
  # set, but Singular.mstd still provides a good heuristic to find a small
  # generating set.
  Q = base_ring(I)
  # Temporary workaround, see #3499
  unique!(filter!(!iszero, Q.(small_generating_set(saturated_ideal(I); algorithm))))
end

# in the graded case, reusing a cached gb makes sense, so use minimal_generating set there
small_generating_set(I::MPolyQuoIdeal{<:MPolyDecRingElem}) = minimal_generating_set(I)

################################################################################
#
#  Promote rule
#
################################################################################

function AbstractAlgebra.promote_rule(::Type{MPolyQuoRingElem{S}}, ::Type{MPolyQuoRingElem{S}}) where {S <: RingElem}
  return MPolyQuoRingElem{S}
end

function AbstractAlgebra.promote_rule(::Type{MPolyQuoRingElem{S}}, ::Type{T}) where {S, T <: RingElem}
  if AbstractAlgebra.promote_rule(S, T) === S
    return MPolyQuoRingElem{S}
  else
    return Union{}
  end
end

@attr Bool function _is_integral_domain(A::MPolyQuoRing)
  return is_prime(modulus(A))
end

# extension of common functionality
symbols(A::MPolyQuoRing) = symbols(base_ring(A))

########################################################################
# Tensor products of rings
########################################################################

function tensor_product(A::MPolyRing, B::MPolyRing; use_product_ordering::Bool = false)
  kk = coefficient_ring(A)
  @assert kk === coefficient_ring(B) "coefficient rings do not coincide"
  res, a, b = polynomial_ring(kk, symbols(A), symbols(B); cached=false)
  if use_product_ordering
    o = induce(gens(res)[1:ngens(A)], default_ordering(A))*induce(gens(res)[ngens(A) + 1:end], default_ordering(B))
    set_default_ordering!(res, o)
  end
  return res, hom(A, res, a), hom(B, res, b)
end

function tensor_product(A::MPolyRing, B::MPolyQuoRing; use_product_ordering::Bool = false)
  R = base_ring(B)
  AR, inc_A, inc_R = tensor_product(A, R)
  I = ideal(AR, inc_R.(gens(modulus(B))))
  res, pr = quo(AR, I)
  if use_product_ordering
    o = induce(gens(AR)[1:ngens(A)], default_ordering(A))*induce(gens(AR)[ngens(A) + 1:end], default_ordering(B))
    set_default_ordering!(AR, o)
  end
  return res, hom(A, res, pr.(inc_A.(gens(A))); check = false), hom(B, res, pr.(inc_R.(gens(R))); check = false)
end

function tensor_product(A::MPolyQuoRing, B::MPolyQuoRing; use_product_ordering::Bool = false)
  RA = base_ring(A)
  RB = base_ring(B)
  P, inc_A, inc_B = tensor_product(RA, RB)
  I = ideal(P, inc_A.(gens(modulus(A)))) + ideal(P, inc_B.(gens(modulus(B))))
  res, pr = quo(P, I)
  if use_product_ordering
    o = induce(gens(P)[1:ngens(A)], default_ordering(A))*induce(gens(P)[ngens(A) + 1:end], default_ordering(B))
    set_default_ordering!(P, o)
  end
  return res, hom(A, res, pr.(inc_A.(gens(RA))); check=false), hom(B, res, pr.(inc_B.(gens(RB))); check=false)
end

function tensor_product(A::MPolyQuoRing, B::MPolyRing; use_product_ordering::Bool = false)
  RA = base_ring(A)
  P, inc_A, inc_B = tensor_product(RA, B)
  I = ideal(P, inc_A.(gens(modulus(A))))
  res, pr = quo(P, I)
  if use_product_ordering
    o = induce(gens(P)[1:ngens(A)], default_ordering(A))*induce(gens(P)[ngens(A) + 1:end], default_ordering(B))
    set_default_ordering!(P, o)
  end
  return res, hom(A, res, pr.(inc_A.(gens(RA))); check=false), hom(B, res, pr.(inc_B.(gens(B))); check = false)
end
