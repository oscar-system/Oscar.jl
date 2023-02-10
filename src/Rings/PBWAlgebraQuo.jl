export PBWAlgQuo, PBWAlgQuoElem

mutable struct PBWAlgQuo{T, S} <: NCRing
  I::PBWAlgIdeal{0, T, S}
end

mutable struct PBWAlgQuoElem{T, S} <: NCRingElem
  parent::PBWAlgQuo{T, S}
  data::PBWAlgElem{T, S}
end

function is_domain_type(::Type{U}) where {U <: PBWAlgQuoElem}
   return false
end

function is_exact_type(a::Type{U}) where {T, U <: PBWAlgQuoElem{T}}
   return is_exact_type(T)
end

elem_type(::PBWAlgQuo{T, S}) where {T, S} = PBWAlgQuoElem{T, S}

elem_type(::Type{PBWAlgQuo{T, S}}) where {T, S} = PBWAlgQuoElem{T, S}

parent_type(::PBWAlgQuoElem{T, S}) where {T, S} = PBWAlgQuo{T, S}

parent_type(::Type{PBWAlgQuoElem{T, S}}) where {T, S} = PBWAlgQuo{T, S}

parent(a::PBWAlgQuoElem) = a.parent

symbols(Q::PBWAlgQuo) = symbols(Q.I.basering)

coefficient_ring(Q::PBWAlgQuo) = coefficient_ring(Q.I.basering)

coefficient_ring(a::PBWAlgQuoElem) = coefficient_ring(parent(a))

modulus(Q::PBWAlgQuo) = Q.I

base_ring(Q::PBWAlgQuo) = Q.I.basering

base_ring(a::PBWAlgQuoElem) = base_ring(parent(a))

function Base.deepcopy_internal(a::PBWAlgQuoElem, dict::IdDict)
  return PBWAlgQuoElem(parent(a), deepcopy_internal(a.data, dict))
end

function expressify(a::PBWAlgQuoElem; context = nothing)
  return expressify(a.data; context=context)
end

@enable_all_show_via_expressify PBWAlgQuoElem

function expressify(Q::PBWAlgQuo; context = nothing)
  return Expr(:call, :/, expressify(Q.I.basering; context = nothing),
                         expressify(Q.I; context = nothing))
end

@enable_all_show_via_expressify PBWAlgQuo

####

function ngens(Q::PBWAlgQuo)
  return ngens(base_ring(Q))
end

function gens(Q::PBWAlgQuo)
  return elem_type(Q)[PBWAlgQuoElem(Q, x) for x in gens(base_ring(Q))]
end

function gen(Q::PBWAlgQuo, i::Int)
  return PBWAlgQuoElem(Q, gen(base_ring(Q), i))
end

function Base.getindex(Q::PBWAlgQuo, i::Int)
  return gen(Q, i)
end

function zero(Q::PBWAlgQuo)
  return PBWAlgQuoElem(Q, zero(base_ring(Q)))
end

function one(Q::PBWAlgQuo)
  return PBWAlgQuoElem(Q, one(base_ring(Q)))
end

function simplify!(a::PBWAlgQuoElem)
  I = parent(a).I
  groebner_assure!(I)
  a.data.sdata = Singular.reduce(a.data.sdata, I.gb)
  return a
end

function Base.hash(a::PBWAlgQuoElem, h::UInt)
  simplify!(a)
  return hash(a.data, h)
end

function iszero(a::PBWAlgQuoElem)
  iszero(a.data) && return true
  I = parent(a).I
  groebner_assure!(I)
  return iszero(Singular.reduce(a.data.sdata, I.gb))
end

function is_unit(a::PBWAlgQuoElem)
  is_unit(a.data) && return true
  is_zero(a.data) && return false
  throw(NotImplementedError(:is_unit, a))
end

function Base.:(==)(a::PBWAlgQuoElem, b::PBWAlgQuoElem)
  return iszero(a - b)
end

function Base.:+(a::PBWAlgQuoElem, b::PBWAlgQuoElem)
  @assert parent(a) == parent(b)
  return PBWAlgQuoElem(parent(a), a.data + b.data)
end

function Base.:-(a::PBWAlgQuoElem, b::PBWAlgQuoElem)
  @assert parent(a) == parent(b)
  return PBWAlgQuoElem(parent(a), a.data - b.data)
end

function Base.:-(a::PBWAlgQuoElem)
  return PBWAlgQuoElem(parent(a), -a.data)
end

function Base.:*(a::PBWAlgQuoElem, b::PBWAlgQuoElem)
  return PBWAlgQuoElem(parent(a), a.data*b.data)
end

function Base.:^(a::PBWAlgQuoElem, b::Int)
  return PBWAlgQuoElem(parent(a), a.data^b)
end

function divexact_left(a::PBWAlgQuoElem, b::PBWAlgQuoElem; check::Bool = true)
  throw(NotImplementedError(:divexact_left, a, b))
end

function divexact_right(a::PBWAlgQuoElem, b::PBWAlgQuoElem; check::Bool = true)
  throw(NotImplementedError(:divexact_right, a, b))
end

####

@doc Markdown.doc"""
    quo(A::PBWAlgRing, I::PBWAlgIdeal)

Given a two-sided ideal `I` of `A`, create the quotient algebra $A/I$ and
return the new algebra together with the quotient map $A\to A/I$.

# Examples
```jldoctest
julia> R, (x, y, z) = QQ["x", "y", "z"];

julia> L = [-x*y, -x*z, -y*z];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x, y, z) = pbw_algebra(R, REL, deglex(gens(R)))
(PBW-algebra over Rational Field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z, PBWAlgElem{fmpq, Singular.n_Q}[x, y, z])

julia> I = two_sided_ideal(A, [x^2, y^2, z^2])
two_sided_ideal(x^2, y^2, z^2)

julia> Q, q = quo(A, I);

julia> Q
(PBW-algebra over Rational Field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z)/two_sided_ideal(x^2, y^2, z^2)

julia> q
Map from
PBW-algebra over Rational Field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z to (PBW-algebra over Rational Field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z)/two_sided_ideal(x^2, y^2, z^2) defined by a julia-function with inverse
```

!!! note
    The example above, shows one way of constructing the exterior algebra on the variables `x`, `y`, `z` over $\mathbb Q$.
    For reasons of efficiency, it is recommended to use the built-in constructor `exterior_algebra` when working with 
    exterior algebras in OSCAR.
"""
function quo(Q::PBWAlgRing, I::PBWAlgIdeal)
  @assert Q == base_ring(I)
  q = PBWAlgQuo(I)
  function im(a::PBWAlgElem)
    @assert parent(a) == Q
    return PBWAlgQuoElem(q, a)
  end
  function pr(a::PBWAlgQuoElem)
    return a.data
  end
  return q, MapFromFunc(im, pr, Q, q)
end

function AbstractAlgebra.promote_rule(::Type{PBWAlgQuoElem{T, S}}, ::Type{PBWAlgQuoElem{T, S}}) where {T, S}
  return PBWAlgQuoElem{T, S}
end

function AbstractAlgebra.promote_rule(::Type{PBWAlgQuoElem{T, S}}, ::Type{U}) where {T, S, U}
  a = AbstractAlgebra.promote_rule(T, U)
  return a == T ? PBWAlgQuoElem{T, S} : Union{}
end

function lift(a::PBWAlgQuoElem)
  return a.data
end

function (Q::PBWAlgQuo)(a::PBWAlgElem)
  return PBWAlgQuoElem(Q, base_ring(Q)(a))
end

function (Q::PBWAlgQuo)()
  return PBWAlgQuoElem(Q, base_ring(Q)())
end

function (Q::PBWAlgQuo{T, S})(c::T) where {T, S}
  return PBWAlgQuoElem(Q, base_ring(Q)(c))
end

function (Q::PBWAlgQuo)(c::IntegerUnion)
  return PBWAlgQuoElem(Q, base_ring(Q)(c))
end

function (Q::PBWAlgQuo)(a::PBWAlgQuoElem)
  parent(a) == Q || error("coercion impossible")
  return a
end

#############################################
##### Exterior algebras
#############################################

@doc Markdown.doc"""
    exterior_algebra(K::Ring, xs::Union{AbstractVector{<:AbstractString}, 
                                    AbstractVector{Symbol}, AbstractVector{Char}})

Given a field `K` and a vector `xs` of,  say, $n$ Strings, Symbols, or Characters, return the $n$-th exterior algebra over `K`.

The generators of the returned algebra print according to the entries of `xs`. See the example below.

# Examples
"""
function exterior_algebra(K::Ring, xs::Union{AbstractVector{<:AbstractString}, 
                                    AbstractVector{Symbol}, AbstractVector{Char}})
  throw(NotImplementedError(:exterior_algebra, K, xs))
end
