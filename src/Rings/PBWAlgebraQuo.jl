export PBWAlgQuo, PBWAlgQuoElem

# DEVELOPER DOC  (code revised 2023-03-early)
#  The revised code allows 2 slightly different implementations:
#  (1) original impl where PBWAlgQuoElem is just a PBWAlgElem with a
#      simple (almost transparent) wrapper -- all arithmetic is
#      delegated to the PBWAlg (hence values are not automatic reduced
#      after multiplication).
#
#  (2) modified impl, currently just for the case of exterior algebra.
#      Arithmetic is delegated to the Singular impl of exterior algebra
#      (which does automatically reduced after multiplication).  To
#      achieve this a "new" pointer to the Singular ring is needed;
#      this pointer is a new data field (sring) in PBQAlgQuo.
#
#  To preserve the original behavior the new field "sdata" in PBWAlgQuo
#  is set to the same value as the field "sdata" in PBWAlg (unless
#  created by constructor for exterior_algebra.
#
# NOTE ABOUT PBWAlgQuoElem
#   Values of type PBWAlgQuoElem are "wrapped" many times: the situation is
#       an internal Singular value: a NC-polynomial & pointer to "quotient ring"
#                                   (which contains the overlying PBWAlg in Singular
#                                    & the reducing ideal)
#       The NC-polynomial has a data repr which is compatible with the overlying
#       Singular PBWAlg, so may be operated upon also by that ring (inside Singular).
#       The Singular repr is then "wrapped" into an Oscar PBWAlgElem (which also
#       contains a "parent" reference to the Oscar object representing the PBWAlg).
#       This PBWAlgElem is "wrapped" once again into  PBWAlgQuoElem (which also
#       contains a "parent" reference to the Oscar object representing the PBWAlgQuo).
#   EXERCISE: sit down and draw a box-and-arrow diagram of values and references!!



###### @attributes    ### DID NOT WORK -- alternative solution found (via has_special_impl, see ExteriorAlgebra.jl)
@attributes mutable struct PBWAlgQuo{T, S} <: NCRing
    I::PBWAlgIdeal{0, T, S}
    sring::Singular.PluralRing{S}  # For ExtAlg this is the Singular impl; o/w same as I.basering.sring

    function PBWAlgQuo(I::PBWAlgIdeal{0, T, S}, sring::Singular.PluralRing{S}) where {T, S}
        return new{T, S}(I, sring)
    end
end


# For backward compatibility: ctor with 1 arg:
#    uses "default" arith impl -- namely that from basering!
function PBWAlgQuo(I::PBWAlgIdeal{0, T, S})  where {T, S}
    return PBWAlgQuo(I, I.basering.sring)
end


function has_special_impl(Q::PBWAlgQuo)
    return (Q.sring != Q.I.basering.sring)
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

elem_type(::Type{PBWAlgQuo{T, S}}) where {T, S} = PBWAlgQuoElem{T, S}

parent_type(::Type{PBWAlgQuoElem{T, S}}) where {T, S} = PBWAlgQuo{T, S}

parent(a::PBWAlgQuoElem) = a.parent

symbols(Q::PBWAlgQuo) = symbols(Q.I.basering)  # EQUIV symbols(Q.sring)  ???

coefficient_ring(Q::PBWAlgQuo) = coefficient_ring(Q.I.basering)  # EQUIV coefficient_ring(Q.sring)  ???

coefficient_ring(a::PBWAlgQuoElem) = coefficient_ring(parent(a))

modulus(Q::PBWAlgQuo) = Q.I

base_ring(Q::PBWAlgQuo) = base_ring(Q.I)

base_ring_type(::Type{PBWAlgQuo{T, S}}) where {T, S} = base_ring_type(PBWAlgIdeal{0, T, S})

function Base.deepcopy_internal(a::PBWAlgQuoElem, dict::IdDict)
  return PBWAlgQuoElem(parent(a), deepcopy_internal(a.data, dict))
end

function expressify(a::PBWAlgQuoElem; context = nothing)
  return expressify(a.data; context=context)
end

@enable_all_show_via_expressify PBWAlgQuoElem

function Base.show(io::IO, Q::PBWAlgQuo)
  @show_name(io, Q)
  @show_special(io, Q)
  print(io, "(", base_ring(Q), ")/", modulus(Q))
end

####

function number_of_generators(Q::PBWAlgQuo)
  return number_of_generators(base_ring(Q))  # EQUIV  number_of_generators(Q.sring)  ???
end

function gens(Q::PBWAlgQuo)
  return elem_type(Q)[PBWAlgQuoElem(Q, PBWAlgElem(Q.I.basering,x)) for x in gens(Q.sring)]
end

function gen(Q::PBWAlgQuo, i::Int)
  return PBWAlgQuoElem(Q, PBWAlgElem(Q.I.basering, gen(Q.sring, i)))
end

function zero(Q::PBWAlgQuo)
  return PBWAlgQuoElem(Q, PBWAlgElem(Q.I.basering, zero(Q.sring)))
end

function one(Q::PBWAlgQuo)
  return PBWAlgQuoElem(Q, PBWAlgElem(Q.I.basering, one(Q.sring)))
end

function simplify(a::PBWAlgQuoElem)
    if has_special_impl(parent(a))
        return a   # short-cut for impls with reducing arithmetic (e.g. exterior algebras)
    end
    I = parent(a).I
    groebner_assure!(I)
    a.data.sdata = Singular.reduce(a.data.sdata, I.gb)
    return a
end

function Base.hash(a::PBWAlgQuoElem, h::UInt)
  simplify(a)
  return hash(a.data, h)
end

function is_zero(a::PBWAlgQuoElem)
    if !has_special_impl(parent(a))  # must reduce if not exterior algebras
        simplify(a)  # see GitHub discussion #2014 -- is_zero can modify repr of its arg!
    end
    return is_zero(a.data.sdata)  # EQUIV  is_zero(a.data)
end

function is_unit(a::PBWAlgQuoElem)
  is_unit(a.data) && return true
  is_zero(a.data) && return false
  throw(NotImplementedError(:is_unit, a))
end

function Base.:(==)(a::PBWAlgQuoElem, b::PBWAlgQuoElem)
  return is_zero(a - b)
end

function Base.:+(a::PBWAlgQuoElem, b::PBWAlgQuoElem)
  check_parent(a, b)
  return PBWAlgQuoElem(parent(a), a.data + b.data)
end

function Base.:-(a::PBWAlgQuoElem, b::PBWAlgQuoElem)
  check_parent(a, b)
  return PBWAlgQuoElem(parent(a), a.data - b.data)
end

function Base.:-(a::PBWAlgQuoElem)
  return PBWAlgQuoElem(parent(a), -a.data)
end

function Base.:*(a::PBWAlgQuoElem, b::PBWAlgQuoElem)
  check_parent(a, b)
  return PBWAlgQuoElem(parent(a), a.data * b.data)
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

@doc raw"""
    quo(A::PBWAlgRing, I::PBWAlgIdeal)

Given a two-sided ideal `I` of `A`, create the quotient algebra $A/I$ and
return the new algebra together with the quotient map $A\to A/I$.

# Examples
```jldoctest
julia> R, (x, y, z) = QQ[:x, :y, :z];

julia> L = [-x*y, -x*z, -y*z];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x, y, z) = pbw_algebra(R, REL, deglex(gens(R)))
(PBW-algebra over rational field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z, PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, z])

julia> I = two_sided_ideal(A, [x^2, y^2, z^2])
two_sided_ideal(x^2, y^2, z^2)

julia> Q, q = quo(A, I);

julia> Q
(PBW-algebra over rational field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z)/two_sided_ideal(x^2, y^2, z^2)

julia> q
Map defined by a julia-function with inverse
  from PBW-algebra over rational field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z
  to (PBW-algebra over rational field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z)/two_sided_ideal(x^2, y^2, z^2)
```

!!! note
    The example above, shows one way of constructing the exterior algebra on the variables `x`, `y`, `z` over $\mathbb Q$.
    For reasons of efficiency, it is recommended to use the built-in constructor `exterior_algebra` when working with 
    exterior algebras in OSCAR.
"""
function quo(Q::PBWAlgRing, I::PBWAlgIdeal; special_impl::Union{Nothing, Singular.PluralRing} = nothing)
  @assert Q == base_ring(I)
  # we assume that special_impl is correct, if given
  if isnothing(special_impl)
    q = PBWAlgQuo(I)
  else
    q = PBWAlgQuo(I, special_impl)
  end
  function im(a::PBWAlgElem)
    @assert parent(a) == Q
    return PBWAlgQuoElem(q, a)
  end
  function pr(a::PBWAlgQuoElem)
    return a.data
  end
  return q, MapFromFunc(Q, q, im, pr)
end


# For some reason we need to specify a promote_rule when both types are the same (why?)
# Doc might be in packages/AbstractAlgebra/*/docs/src/ring_interface.md (near line 580)
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
  return PBWAlgQuoElem(Q, PBWAlgElem(base_ring(Q), Q.sring(0)))
end

function (Q::PBWAlgQuo{T, S})(c::T) where {T, S}
  return PBWAlgQuoElem(Q, PBWAlgElem(base_ring(Q), Q.sring(c)))
end

function (Q::PBWAlgQuo)(c::IntegerUnion)
  return PBWAlgQuoElem(Q, PBWAlgElem(base_ring(Q), Q.sring(c)))
end

function (Q::PBWAlgQuo)(a::PBWAlgQuoElem)
  @req parent(a) == Q "coercion between different PBWAlg quotients not possible"
  return a
end

### Conformance test element generation
function ConformanceTests.generate_element(Q::PBWAlgQuo{QQFieldElem})
  R = base_ring(Q)
  return Q(R(rand(base_ring(R), 1:4, 1:4, 1:4)))
end

#############################################
##### Exterior algebras
#############################################

#---------------------------------------------------------
# SEE FILE experimental/ExteriorAlgebra/ExteriorAlgebra.jl
#---------------------------------------------------------


# 2023-03-09 JAA  commented out placeholder code below -- should be replaced by code from ExteriorAlgebra.jl

# @doc raw"""
#     exterior_algebra(K::Ring, xs::AbstractVector{<:VarName})

# Given a field `K` and a vector `xs` of,  say, $n$ Strings, Symbols, or Characters, return the $n$-th exterior algebra over `K`.

# The generators of the returned algebra print according to the entries of `xs`. See the example below.

# # Examples
# """
# function exterior_algebra(K::Ring, xs::AbstractVector{<:VarName})
#   throw(NotImplementedError(:exterior_algebra, K, xs))
# end
