###############################################################################
#
#  Abelian closure of the rationals
#
###############################################################################

# This is an implementation of Q^ab, the abelian closure of the rationals,
# which is modeled as the union of cyclotomic fields.
#
# We make Q^ab a singleton, similar to ZZ and QQ. Thus there will ever be only
# one copy of Q^ab. In particular the elements do not have a parent stored.
# One reason for this decision is that in the circumstances where we decide
# to construct elements of Q^ab, we do not want the user to supply this field.
# This is the case for character tables as well as functions related to
# binomial ideals.
#
# This has the minor problem that printing is controlled using global state.
#
# Note that there are two possibilities construct a nth root of unity when n is
# even and n%4!=0. either we can construct the field Q(z_n) or we take -z_(n/2)
# as a primitive n-th root. to change between these two options, use
# PCharSaturateAll with allroots or allrootsNew (change this in the code)

abstract type CyclotomicField end

export CyclotomicField

module AbelianClosure 

using ..Oscar

import Base: +, *, -, //, ==, zero, one, ^, div, isone, iszero,
             deepcopy_internal, hash, reduce

#import ..Oscar.AbstractAlgebra: promote_rule

import ..Oscar: AbstractAlgebra, addeq!, characteristic, elem_type, divexact, gen,
                has_preimage_with_preimage, is_root_of_unity, is_unit, mul!, parent,
                parent_type, promote_rule, root, root_of_unity, roots

using Hecke
import Hecke: conductor, data

################################################################################
#
#  Types
#
################################################################################

@attributes mutable struct QQAbField{T} <: Nemo.Field # union of cyclotomic fields
  fields::Dict{Int, T} # Cache for the cyclotomic fields
  s::String

  function QQAbField{T}(fields::Dict{Int, T}) where T
    return new(fields)
  end
end

const _QQAb = QQAbField{AbsSimpleNumField}(Dict{Int, AbsSimpleNumField}())
const _QQAb_sparse = QQAbField{AbsNonSimpleNumField}(Dict{Int, AbsNonSimpleNumField}())

mutable struct QQAbElem{T} <: Nemo.FieldElem
  data::T                             # Element in cyclotomic field
  c::Int                              # Conductor of field
end
#T test that data really belongs to a cyclotomic field!

# This is a functor like object G with G(n) = primitive n-th root of unity

mutable struct QQAbFieldGen{T}
  K::QQAbField{T}
end

const _QQAbGen = QQAbFieldGen(_QQAb)
const _QQAbGen_sparse = QQAbFieldGen(_QQAb_sparse)

################################################################################
#
#  Creation of the field
#
################################################################################

@doc raw"""
    abelian_closure(QQ::RationalField)

Return a pair `(K, z)` consisting of the abelian closure `K` of the rationals
and a generator `z` that can be used to construct primitive roots of unity in
`K`.

An optional keyword argument `sparse` can be set to `true` to switch to a 
sparse representation. Depending on the application this can be much faster
or slower.

# Examples
```jldoctest; setup = :(using Oscar)
julia> K, z = abelian_closure(QQ);

julia> z(36)
zeta(36)

julia> K, z = abelian_closure(QQ, sparse = true);

julia> z(36)
-zeta(36, 9)*zeta(36, 4)^4 - zeta(36, 9)*zeta(36, 4)

```
"""
function abelian_closure(::QQField; sparse::Bool = false) 
  if sparse 
    return _QQAb_sparse, _QQAbGen_sparse
  else
    return _QQAb, _QQAbGen
  end
end

"""
    gen(K::QQAbField)

Return the generator of the abelian closure `K` that can be used to construct
primitive roots of unity.
"""
gen(K::QQAbField{AbsSimpleNumField}) = _QQAbGen
gen(K::QQAbField{AbsNonSimpleNumField}) = _QQAbGen_sparse

"""
    gen(K::QQAbField, s::String)

Return the generator of the abelian closure `K` that can be used to construct
primitive roots of unity. The string `s` will be used during printing.
"""
function gen(K::QQAbField, s::String)
  K.s = s
  return gen(K)
end

function characteristic(::QQAbField)
  return 0
end

################################################################################
#
#  Parent and element functions
#
################################################################################

elem_type(::Type{QQAbField{AbsSimpleNumField}}) = QQAbElem{AbsSimpleNumFieldElem}
parent_type(::Type{QQAbElem{AbsSimpleNumFieldElem}}) = QQAbField{AbsSimpleNumField}
parent(::QQAbElem{AbsSimpleNumFieldElem}) = _QQAb

elem_type(::Type{QQAbField{AbsNonSimpleNumField}}) = QQAbElem{AbsNonSimpleNumFieldElem}
parent_type(::Type{QQAbElem{AbsNonSimpleNumFieldElem}}) = QQAbField{AbsNonSimpleNumField}
parent(::QQAbElem{AbsNonSimpleNumFieldElem}) = _QQAb_sparse

################################################################################
#
#  Field access
#
################################################################################

function _variable(K::QQAbField)
  if isdefined(K, :s)
    return K.s
  elseif Oscar.is_unicode_allowed()
    return "ζ"
  else
    return "zeta"
  end
end

_variable(b::QQAbElem{AbsSimpleNumFieldElem}) = Expr(:call, Symbol(_variable(_QQAb)), b.c)

function _variable(b::QQAbElem{AbsNonSimpleNumFieldElem}) 
  k = parent(b.data)
  lc = get_attribute(k, :decom)
  n = get_attribute(k, :cyclo)
  return [Expr(:call, Symbol(_variable(parent(b))), n, divexact(n, i)) for i = lc]
end

function Hecke.cyclotomic_field(K::QQAbField{AbsSimpleNumField}, c::Int)
  if haskey(K.fields, c)
    k = K.fields[c]
    return k, gen(k)
  else
    k, z = cyclotomic_field(c, string("\$", "(", c, ")"), cached = false)
    K.fields[c] = k
    return k, z
  end
end

function ns_gen(K::AbsNonSimpleNumField)
  #z_pq^p = z_q and z_pg^q = z_p
  #thus z_pq = z_p^a z_q^b implies
  #z_pq^p = z_q^pb, so pb = 1 mod q
  #so:
  lc = get_attribute(K, :decom)
  n = get_attribute(K, :cyclo)
  return prod(gen(K, i)^invmod(divexact(n, lc[i]), lc[i]) for i=1:length(lc))
end

function Hecke.cyclotomic_field(K::QQAbField{AbsNonSimpleNumField}, c::Int)
  if haskey(K.fields, c)
    k = K.fields[c]
    return k, ns_gen(k)
  else
    k, _ = cyclotomic_field(NonSimpleNumField, c, "\$")
    K.fields[c] = k
    return k, ns_gen(k)
  end
end

Hecke.data(a::QQAbElem) = a.data

################################################################################
#
#  Creation of elements
#
################################################################################

# This function finds a primitive root of unity in our field, note this is
# not always e^(2*pi*i)/n

function root_of_unity(K::QQAbField{AbsSimpleNumField}, n::Int)
  if n % 2 == 0 && n % 4 != 0
    c = div(n, 2)
  else
    c = n
  end
  K, z = cyclotomic_field(K, c)
  if c == n
    return QQAbElem{AbsSimpleNumFieldElem}(z, c)
  else
    return QQAbElem{AbsSimpleNumFieldElem}(-z, c)
  end
end

function root_of_unity(K::QQAbField{AbsNonSimpleNumField}, n::Int)
  if n % 2 == 0 && n % 4 != 0
    c = div(n, 2)
  else
    c = n
  end
  K, z = cyclotomic_field(K, c)
  if c == n
    return QQAbElem{AbsNonSimpleNumFieldElem}(z, c)
  else
    return QQAbElem{AbsNonSimpleNumFieldElem}(-z, c)
  end
end


function root_of_unity2(K::QQAbField, c::Int)
  # This function returns the primitive root of unity e^(2*pi*i/n)
  K, z = cyclotomic_field(K, c)
  return QQAbElem(z, c)
end

(z::QQAbFieldGen)(n::Int) = root_of_unity(z.K, n)
(z::QQAbFieldGen)(n::Int, r::Int) = z(n)^r

one(K::QQAbField) = K(1)

one(a::QQAbElem) = one(parent(a))

function isone(a::QQAbElem)
  return isone(data(a))
end

function iszero(a::QQAbElem)
  return iszero(data(a))
end

zero(K::QQAbField) = K(0)

zero(a::QQAbElem) = zero(parent(a))

function (K::QQAbField)(a::Union{ZZRingElem, QQFieldElem, Integer, Rational})
  return a*root_of_unity(K, 1)
end

function (K::QQAbField)(a::QQAbElem)
  return a
end

(K::QQAbField)() = zero(K)

function (K::QQAbField)(a::AbsSimpleNumFieldElem)
  F = parent(a)

  # Cyclotomic fields are naturally embedded into `K`.
  fl, f = Hecke.is_cyclotomic_type(F)
  fl && return QQAbElem(a, f)

  # Quadratic fields are naturally embedded into `K`.
  if degree(F) == 2
    # If the defining polynomial of `F` is `X^2 + A X + B` then
    # `D = A^2 - 4 B` is a square in `F` (cf. [Coh93, p. 218]).
    pol = F.pol
    A = coeff(pol, 1)
    D = A^2 - 4*coeff(pol, 0)
    Dn = numerator(D)
    Dd = denominator(D)
    f = Dn * Dd

    x = coeff(a, 0)
    y = coeff(a, 1)
    iszero(y) && return QQAbElem(parent(a)(x), 1)
    d = sign(f)
    c = 1
    for (p, e) in factor(f)
      if mod(e, 2) == 1
        d = d*p
        c = c*p^div(e-1, 2)
      else
        c = c*p^div(e, 2)
      end
    end
    # `d` is the signed squarefree part of `f`, and `f == d*c^2`
    if mod(d, 4) == 1
      N = abs(d)
    else
      N = 4*abs(d)
    end
    r = square_root_in_cyclotomic_field(K, Int(d), Int(N))
    return (x - y*A//2) + (y*c//(2*Dd)) * r
  end

  # We have no natural embeddings for other (abelian) number fields.
  throw(ArgumentError("no natural embedding of $(parent(a)) into QQAbField"))
end

################################################################################
#
#  String I/O
#
################################################################################

function Base.show(io::IO, a::QQAbField{AbsNonSimpleNumField})
  print(io, "(Sparse) abelian closure of Q")
end

function Base.show(io::IO, a::QQAbField{AbsSimpleNumField})
  print(io, "Abelian closure of Q")
end

function Base.show(io::IO, a::QQAbFieldGen)
  if isa(a.K, QQAbField{AbsSimpleNumField})
    print(io, "Generator of abelian closure of Q")
  else
    print(io, "Generator of sparse abelian closure of Q")
  end
end

"""
    set_variable!(K::QQAbField, s::String)

Change the printing of the primitive n-th root of the abelian closure of the
rationals to `s(n)`, where `s` is the supplied string.
"""
function set_variable!(K::QQAbField, s::String)
  ss = _variable(K)
  K.s = s
  return ss
end

"""
    get_variable(K::QQAbField)

Return the string used to print the primitive n-th root of the abelian closure
of the rationals.
"""
get_variable(K::QQAbField) = _variable(K)

function AbstractAlgebra.expressify(b::QQAbElem{AbsSimpleNumFieldElem}; context = nothing)
  a = data(b)
  return AbstractAlgebra.expressify(parent(parent(a).pol)(a), _variable(b), context = context)
end

function AbstractAlgebra.expressify(b::QQAbElem{AbsNonSimpleNumFieldElem}; context = nothing)
  a = data(b)
  return AbstractAlgebra.expressify(a.data, _variable(b), context = context)
end

Oscar.@enable_all_show_via_expressify QQAbElem

################################################################################
#
#  Singular ring
#
################################################################################

function Oscar.singular_coeff_ring(F::QQAbField)
  return Singular.CoefficientRing(F)
end

################################################################################
#
#  Coercion between cyclotomic fields
#
################################################################################

function is_conductor(n::Int)
  if isodd(n)
    return true
  end
  return n % 4 == 0
end

function coerce_up(K::AbsSimpleNumField, n::Int, a::QQAbElem{AbsSimpleNumFieldElem})
  d = div(n, a.c)
  @assert n % a.c == 0
  #z_n^(d) = z_a
  R = parent(parent(data(a)).pol)
  return QQAbElem{AbsSimpleNumFieldElem}(evaluate(R(data(a)), gen(K)^d), n)
end

function coerce_up(K::AbsNonSimpleNumField, n::Int, a::QQAbElem{AbsNonSimpleNumFieldElem})
  d = div(n, a.c)
  @assert n % a.c == 0
  lk = get_attribute(parent(a.data), :decom)
  #gen(k, i) = gen(K, j)^n for the unique j s.th. gcd(lk[i], lK[j])
  # and n = lK[j]/lk[i]
  #z_n^(d) = z_a
  return QQAbElem{AbsNonSimpleNumFieldElem}(evaluate(data(a).data, [ns_gen(K)^divexact(n, i) for i=lk]), n)
end


function coerce_down(K::AbsSimpleNumField, n::Int, a::QQAbElem)
  throw(Hecke.NotImplemented())
end

function make_compatible(a::QQAbElem{T}, b::QQAbElem{T}) where {T}
  if a.c == b.c
    return a,b
  end
  d = lcm(a.c, b.c)
  K, = cyclotomic_field(parent(a), d)
  return coerce_up(K, d, a), coerce_up(K, d, b)
end

function minimize(::typeof(CyclotomicField), a::AbstractArray{AbsSimpleNumFieldElem})
  fl, c = Hecke.is_cyclotomic_type(parent(a[1]))
  @assert all(x->parent(x) == parent(a[1]), a)
  @assert fl
  for p = keys(factor(c).fac)
    while c % p == 0
      K, _ = cyclotomic_field(Int(div(c, p)), cached = false)
      b = similar(a)
      OK = true
      for x = eachindex(a)
        y = Hecke.force_coerce_cyclo(K, a[x], Val{false})
        if y === nothing
          OK = false
        else
          b[x] = y
        end
      end
      if OK
        a = b
        c = div(c, p)
      else
        break
      end
    end
  end
  return a
end

function minimize(::typeof(CyclotomicField), a::MatElem{AbsSimpleNumFieldElem})
  return matrix(minimize(CyclotomicField, a.entries))
end

function minimize(::typeof(CyclotomicField), a::AbsSimpleNumFieldElem)
  return minimize(CyclotomicField, [a])[1]
end

#TODO:
# Here we use conductor in the sense that
# an abelian number field K has conductor n iff the n-th cyclotomic field
# is the smallest cyclotomic field that contains K,
# and the conductor of a field element is the conductor of the field
# it generates.
# Claus says that the conductor of a field element can also be read
# w.r.t. an order.
# Do we have a naming problem?
# (If not then we can just add documentation.)
conductor(a::AbsSimpleNumFieldElem) = conductor(parent(minimize(CyclotomicField, a)))

function conductor(k::AbsSimpleNumField)
  f, c = Hecke.is_cyclotomic_type(k)
  f || error("field is not of cyclotomic type")
  return c
end

conductor(a::QQAbElem) = conductor(data(a))

# What we want is the conductor of the domain of the map, but we need the map.
function conductor(phi::MapFromFunc{T, QQAbField{T}}) where T
  return lcm([conductor(phi(x)) for x in gens(domain(phi))])
end

################################################################################
#
#  Conversions to `ZZRingElem` and `QQFieldElem` (like for `AbsSimpleNumFieldElem`)
#
################################################################################

(R::QQField)(a::QQAbElem) = R(a.data)
(R::ZZRing)(a::QQAbElem) = R(a.data)


################################################################################
#
#  Ring interface functions
#
################################################################################

is_unit(a::QQAbElem) = !iszero(a)

canonical_unit(a::QQAbElem) = a

################################################################################
#
#  Minimal polynomial
#
################################################################################

Hecke.minpoly(a::QQAbElem) = minpoly(data(a))

################################################################################
#
#  Syntactic sugar
#
################################################################################

function Hecke.number_field(::QQField, a::QQAbElem; cached::Bool = false)
  f = minpoly(a)
  k, b = number_field(f, check = false, cached = cached)
  return k, b
end

function Hecke.number_field(::QQField, a::AbstractVector{<: QQAbElem}; cached::Bool = false)
  if length(a) == 0
    return Hecke.rationals_as_number_field()[1]
  end
  f = lcm([Hecke.is_cyclotomic_type(parent(data(x)))[2] for x = a])
  K = cyclotomic_field(f)[1]
  k, mkK = Hecke.subfield(K, [K(data(x)) for x = a])
  return k, gen(k)
end

Base.getindex(::QQField, a::QQAbElem) = number_field(QQ, a)
Base.getindex(::QQField, a::Vector{QQAbElem{T}}) where T = number_field(QQ, a)
Base.getindex(::QQField, a::QQAbElem...) = number_field(QQ, [x for x in a])

################################################################################
#
#  Arithmetic
#
################################################################################

function +(a::QQAbElem, b::QQAbElem)
  a, b = make_compatible(a, b)
  return QQAbElem(data(a) + data(b), a.c)
end

function *(a::QQAbElem, b::QQAbElem)
  a, b = make_compatible(a, b)
  return QQAbElem(data(a) * data(b), a.c)
end

function -(a::QQAbElem)
  return QQAbElem(-data(a), a.c)
end

function ^(a::QQAbElem, n::Integer)
  return QQAbElem(data(a)^n, a.c)
end

function ^(a::QQAbElem, n::ZZRingElem)
  return a^Int(n)
end

function -(a::QQAbElem, b::QQAbElem)
  a, b = make_compatible(a, b)
  return QQAbElem(a.data-b.data, a.c)
end

function //(a::QQAbElem, b::QQAbElem)
  a, b = make_compatible(a, b)
  return QQAbElem(a.data//b.data, a.c)
end

function div(a::QQAbElem, b::QQAbElem)
  a, b = make_compatible(a, b)
  return QQAbElem(a.data//b.data, a.c)
end

function divexact(a::QQAbElem, b::QQAbElem; check::Bool = true)
  a, b = make_compatible(a, b)
  return QQAbElem(divexact(a.data, b.data), a.c)
end

function inv(a::QQAbElem)
  return QQAbElem(inv(data(a)), a.c)
end

################################################################################
#
#  Unsafe operations
#
################################################################################

function addeq!(c::QQAbElem, a::QQAbElem)
  _c, _a = make_compatible(c, a)
  addeq!(_c.data, _a.data)
  return _c
end

function neg!(a::QQAbElem)
  mul!(a.data,a.data,-1)
  return a
end

function mul!(c::QQAbElem, a::QQAbElem, b::QQAbElem)
  a, b = make_compatible(a, b)
  b, c = make_compatible(b, c)
  a, b = make_compatible(a, b)
  mul!(c.data, a.data, b.data)
  return c
end

################################################################################
#
#  Ad hoc binary operations
#
################################################################################

*(a::ZZRingElem, b::QQAbElem) = QQAbElem(b.data*a, b.c)

*(a::QQFieldElem, b::QQAbElem) = QQAbElem(b.data*a, b.c)

*(a::Integer, b::QQAbElem) = QQAbElem(data(b) * a, b.c)

*(a::Rational, b::QQAbElem) = QQAbElem(data(b) * a, b.c)

*(a::QQAbElem, b::ZZRingElem) = b*a

*(a::QQAbElem, b::QQFieldElem) = b*a

*(a::QQAbElem, b::Integer) = b*a

*(a::QQAbElem, b::Rational) = b*a

+(a::ZZRingElem, b::QQAbElem) = QQAbElem(b.data + a, b.c)

+(a::QQFieldElem, b::QQAbElem) = QQAbElem(b.data + a, b.c)

+(a::Integer, b::QQAbElem) = QQAbElem(data(b) + a, b.c)

+(a::Rational, b::QQAbElem) = QQAbElem(data(b) + a, b.c)

+(a::QQAbElem, b::ZZRingElem) = b + a

+(a::QQAbElem, b::QQFieldElem) = b + a

+(a::QQAbElem, b::Integer) = b + a

+(a::QQAbElem, b::Rational) = b + a

-(a::ZZRingElem, b::QQAbElem) = QQAbElem(-(a, data(b)), b.c)

-(a::QQFieldElem, b::QQAbElem) = QQAbElem(-(a, data(b)), b.c)

-(a::Integer, b::QQAbElem) = QQAbElem(-(a, data(b)), b.c)

-(a::Rational, b::QQAbElem) = QQAbElem(-(a, data(b)), b.c)

-(a::QQAbElem, b::ZZRingElem) = QQAbElem(-(data(a), b), a.c)

-(a::QQAbElem, b::QQFieldElem) = QQAbElem(-(data(a), b), a.c)

-(a::QQAbElem, b::Integer) = QQAbElem(-(data(a), b), a.c)

-(a::QQAbElem, b::Rational) = QQAbElem(-(data(a), b), a.c)

//(a::QQAbElem, b::ZZRingElem) = QQAbElem(data(a)//b, a.c)

//(a::QQAbElem, b::QQFieldElem) = QQAbElem(data(a)//b, a.c)

//(a::QQAbElem, b::Integer) = QQAbElem(data(a)//b, a.c)

//(a::QQAbElem, b::Rational) = QQAbElem(data(a)//b, a.c)

################################################################################
#
#  Comparison
#
################################################################################

function ==(a::QQAbElem, b::QQAbElem)
  a, b = make_compatible(a, b)
  return a.data == b.data
end

function ==(a::QQAbElem, b::Union{ZZRingElem, QQFieldElem, Integer, Rational})
  return data(a) == b
end

function ==(a::Union{ZZRingElem, QQFieldElem, Integer, Rational}, b::QQAbElem)
  return b == a
end

hash(a::QQAbElem, h::UInt) = hash(minimize(CyclotomicField, a.data), h)

################################################################################
#
#  Copy and deepcopy
#
################################################################################

function Base.copy(a::QQAbElem)
  return QQAbElem(data(a), a.c)
end

function Base.deepcopy_internal(a::QQAbElem, dict::IdDict)
  return QQAbElem(deepcopy_internal(data(a), dict), a.c)
end

################################################################################
#
#  Promotion rules
#
################################################################################

AbstractAlgebra.promote_rule(::Type{QQAbElem}, ::Type{Int}) = QQAbElem

AbstractAlgebra.promote_rule(::Type{QQAbElem}, ::Type{ZZRingElem}) = QQAbElem

AbstractAlgebra.promote_rule(::Type{QQAbElem}, ::Type{QQFieldElem}) = QQAbElem

###############################################################################
#
#   Functions for computing roots
#
###############################################################################

function Oscar.root(a::QQAbElem, n::Int)
  Hecke.@req is_root_of_unity(a) "Element must be a root of unity"
  o = Oscar.order(a)
  l = o*n
  mu = root_of_unity2(parent(a), Int(l))
  return mu
end

function Oscar.roots(f::PolyRingElem{QQAbElem{T}}) where T
  QQAb = base_ring(f)
  c = reduce(lcm, map(conductor, AbstractAlgebra.coefficients(f)), init = Int(1))
  k, z = cyclotomic_field(QQAb, c)
  f = map_coefficients(x->k(x.data), f)
  lf = factor(f).fac
  #we need to find the correct cyclotomic field...
  #can't use ray_class_group in k as this is expensive (needs class group)
  #need absolute norm group
  QQ = rationals_as_number_field()[1]

  C = cyclotomic_field(ClassField, c)

  rts = QQAbElem{T}[]

  for g = keys(lf)
    c = reduce(lcm, map(conductor, AbstractAlgebra.coefficients(g)), init = Int(1))
    #so THIS factor lives in cyclo(c)
    k, z = cyclotomic_field(QQAb, c)
    d = numerator(norm(k(discriminant(g))))

    R, mR = ray_class_group(lcm(d, c)*maximal_order(QQ), infinite_places(QQ), 
                                        n_quo = degree(g)*degree(k))
    q, mq = quo(R, [R[0]], false)
    for p = PrimesSet(100, -1, c, 1) #totally split primes.
      if d % p == 0
        continue
      end
      if order(q) <= degree(g)*degree(k)
        break
      end
      P = preimage(mR, p*maximal_order(QQ))
      if iszero(mq(P))
        continue
      end

      me = modular_init(k, p)
      lp = Hecke.modular_proj(g, me)
      for pg = lp
        l = factor(pg)
        q, mqq = quo(q, [degree(x)*mq(P) for x = keys(l.fac)], false)
        mq = mq*mqq
        if order(q) <= degree(g)*degree(k)
          break
        end
      end
      if order(q) <= degree(g)*degree(k)
        break
      end
    end
    D = C*ray_class_field(mR, mq)
    c = norm(conductor(D)[1])
    k, a = cyclotomic_field(QQAb, Int(c))
    rt = roots(map_coefficients(k, g))
    append!(rts, map(QQAb, rt))
  end
  return rts::Vector{QQAbElem{T}}
end

function Oscar.roots(a::QQAbElem{T}, n::Int) where {T}
  #strategy:
  # - if a is a root-of-1: trivial, as the answer is also roots-of-1
  # - if a can "easily" be made into a root-of-one: doit
  #   easily is "defined" as <a> = b^n and gens(inv(b))[2]^n*a is a root 
  #   as ideal roots are easy
  # - else: call the function above which is non-trivial...

  corr = one(parent(a))

  if !is_root_of_unity(a) 
    zk = maximal_order(parent(a.data)) #should be for free
    fl, i = is_power(a.data*zk, n)
    _, x = polynomial_ring(parent(a), cached = false)
    fl || return roots(x^n-a)::Vector{QQAbElem{T}}
    b = gens(Hecke.inv(i))[end]
    c = deepcopy(a)
    c.data = b
    corr = Hecke.inv(c)
    a *= c^n
    fl = is_root_of_unity(a)
    fl || return (corr .* roots(x^n-a))::Vector{QQAbElem{T}}
  end
  
  o = order(a)
  l = o*n
  mu = root_of_unity(parent(a), Int(l))
  A = QQAbElem[]
  if l==1 && mu==a
    push!(A, mu)
  end
  for k = 0:(l-1)
    el = mu^k
    if el^n == a
      push!(A, el)
    end
  end
  return [x*corr for x = A]::Vector{QQAbElem{T}}
end

function is_root_of_unity(a::QQAbElem)
  return is_torsion_unit(a.data, true)
  #=
  b = a^a.c
  return b.data == 1 || b.data == -1
  =#
end

function Oscar.order(a::QQAbElem)
  f = Nemo.factor(ZZRingElem(2*a.c))
  o = 1
  for (p, e) = f.fac
    b = a^div(2*a.c, Int(p)^e)
    f = 0

    while !isone(b)
      b = b^p
      f += 1
    end
    o *= p^f
  end
  return o
end

# Convenient sqrt and cbrt functions as simple wrappers around the roots function,
# which is already implemented for QQAbElem directly

function Oscar.sqrt(a::QQAbElem)
  sqrt = Oscar.roots(a, 2)
  if is_empty(sqrt)
    error("Element $a does not have a square root")
  end
  return sqrt[1]
end

function Oscar.cbrt(a::QQAbElem)
  cbrt = Oscar.roots(a,3)
  if is_empty(cbrt)
    error("Element $a does not have a cube root")
  end
  return cbrt[1]
end


###############################################################################
#
#   Embeddings of subfields of cyclotomic fields
#   (works for proper subfields of cycl. fields only if these fields
#   have been constructed as such)
#

# Construct the map from `F` to an abelian closure `K` such that `gen(F)`
# is mapped to `x`.
# If `F` is a cyclotomic field with conductor `N` then assume that `gen(F)`
# is mapped to `QQAbElem(gen(F), N)`.
# (Use that the powers of this element form a basis of the field.)
function _embedding(F::QQField, K::QQAbField{AbsSimpleNumField},
                    x::QQAbElem{AbsSimpleNumFieldElem})
  C1, z = cyclotomic_field(1)

  f = function(x::QQFieldElem)
    return QQAbElem(C1(x), 1)
  end

  finv = function(x::QQAbElem; check::Bool = false)
    if conductor(x) == 1
      return Hecke.force_coerce_cyclo(C1, data(x))
    elseif check
      return
    else
      error("element has no preimage")
    end
  end

  return MapFromFunc(F, K, f, finv)
end

function _embedding(F::AbsSimpleNumField, K::QQAbField{AbsSimpleNumField},
                    x::QQAbElem{AbsSimpleNumFieldElem})
  fl, n = Hecke.is_cyclotomic_type(F)
  if fl
    # This is cheaper.
    f = function(x::AbsSimpleNumFieldElem)
      return QQAbElem(x, n)
    end

    finv = function(x::QQAbElem; check::Bool = false)
      if n % conductor(x) == 0
        return Hecke.force_coerce_cyclo(F, data(x))
      elseif check
        return
      else
        error("element has no preimage")
      end
    end
  else
    # `F` is expected to be a proper subfield of a cyclotomic field.
    n = conductor(x)
    x = data(x)
    Kn, = AbelianClosure.cyclotomic_field(K, n)
    powers = [Hecke.coefficients(Hecke.force_coerce_cyclo(Kn, x^i))
              for i in 0:degree(F)-1]
    c = transpose(matrix(QQ, powers))
    R = parent(F.pol)

    f = function(z::AbsSimpleNumFieldElem)
      return QQAbElem(evaluate(R(z), x), n)
    end

    finv = function(x::QQAbElem; check::Bool = false)
      n % conductor(x) == 0 || return false, zero(F)
      # Write `x` w.r.t. the n-th cyclotomic field ...
      g = gcd(x.c, n)
      Kg, = AbelianClosure.cyclotomic_field(K, g)
      x = Hecke.force_coerce_cyclo(Kg, data(x))
      x = Hecke.force_coerce_cyclo(Kn, x)
      # ... and then w.r.t. `F`
      a = Hecke.coefficients(x)
      fl, sol = can_solve_with_solution(c, matrix(QQ, length(a), 1, a); side = :right)
      if fl
        b = transpose(sol)
        b = [b[i] for i in 1:length(b)]
        return F(b)
      elseif check
        return
      else
        error("element has no preimage")
      end
    end
  end
  return MapFromFunc(F, K, f, finv)
end

# The following works only if `mp.g` admits a second argument,
# which is the case if `mp` has been constructed by `_embedding` above.
function has_preimage_with_preimage(mp::MapFromFunc{AbsSimpleNumField, QQAbField{AbsSimpleNumField}}, x::QQAbElem{AbsSimpleNumFieldElem})
  pre = mp.g(x, check = true)
  if isnothing(pre)
    return false, zero(domain(mp))
  else
    return true, pre
  end
end

###############################################################################
#
#   Galois automorphisms of QQAb
#
################################################################################

# The Galois automorphisms of the $n$-th cyclotomic field are the maps
# defined by $\zeta_n \mapsto \zeta_n^k$, for $1 \leq k < n$,
# with $\gcd(n, k) = 1$.
# Thus we can define automorphisms $\sigma_k$ of QQAb as follows.
# For each prime power $q$, $\zeta_q$ is mapped to $\zeta_q^k$ if
# $k$ and $q$ are coprime, and to $\zeta_q$ otherwise.
#
# The action of such a map $\sigma_k$ on the $n$-th cyclotomic field can be
# described by $\sigma_l$, with $l$ coprime to $n$:
# Write $n = n_0 n_1$ where $\gcd(n_0, n_1) = 1$ and $n_1$ is maximal
# with $\gcd(k, n_1) = 1$, and choose $a, b$ with $1 = a n_0 + b n_1$.
# Then $l = k a n_0 + b n_1$ is coprime to $n$ and has the properties
# $l \equiv 1 \pmod{n_0}$ and $l \equiv k \pmod{n_1}$.

mutable struct QQAbAutomorphism
  exp::Int
end

Oscar.hom(K::QQAbField, L::QQAbField, k::Int) = QQAbAutomorphism(k)

(f::QQAbAutomorphism)(a::QQAbElem) = a^f

function ^(val::QQAbElem, sigma::QQAbAutomorphism)
  k = sigma.exp
  n = val.c
  g = gcd(k, n)
  if g != 1
    # Replace `k` by an equivalent one that is coprime to `n`.
    n0 = 1
    n1 = n
    for (p, exp) in collect(Oscar.factor(g))
      while mod(n1, p) == 0
        n0 = n0*p
        n1 = div(n1, p )
      end
    end
    (gg, a, b) = gcdx(n0, n1)
    @assert gg == 1 "n0 and n1 should be coprime"
    k = k*a*n0 + b*n1
  end
  data = val.data  # AbsSimpleNumFieldElem
  coeffs = Nemo.coefficients(data)
  res = zeros(eltype(coeffs), n)
  res[1] = coeffs[1]
  for i in 2:length(coeffs)
    res[Int(mod((i-1)*k, n)+1)] = coeffs[i]
  end
  F = parent(data) # cycl. field
  R = parent(F.pol)
  return QQAbElem(F(R(res)), n)
end

Base.conj(elm::QQAbElem) = elm^QQAbAutomorphism(-1)

###############################################################################
#
#   Elements in quadratic subfields of cyclotomic fields
#
###############################################################################

function generators_galois_group_cyclotomic_field(n::Int)
  res = GAP.Globals.GeneratorsPrimeResidues(GAP.Obj(n))
  return [QQAbAutomorphism(k)
          for k in Vector{Int}(GAP.Globals.Flat(res.generators))]
end

"""
    square_root_in_cyclotomic_field(F::QQAbField, n::Int, N::Int)

Return an element `a` in the field `F` that is represented w.r.t. the `N`-th
cyclotomic field and has the property `a^2 == n`.

If the `N`-th cyclotomic field does not contain such an element
then `nothing` is returned.

If `n` is positive then `a` is the positive square root of `n`,
otherwise `a` is a positive multiple of the imaginary unit.
(Here we assume that the underlying primitive `N`-th root of unity
is identified with the complex number `exp(2*Pi*i/N)`,
where `i` is the imaginary unit.)
"""
function square_root_in_cyclotomic_field(F::QQAbField, n::Int, N::Int)
  @req N > 0 "conductor ($N) must be positive"
  z = root_of_unity(F, N)
  if n == 0
    return zero(z)
  end

  cf = 1
  sqf = 1
  for (p,e) in collect(factor(n))
    cf = cf * p^div(e, 2)
    if e % 2 != 0
      sqf = sqf * p
    end
  end
  nn = Int(sqf)  # nn is positive and squarefree
  N % nn == 0 || return

  n4 = nn % 4
  if n4 == 1
    if n < 0
      if N % 4 != 0
        return
      end
      N4 = div(N, 4)
      z4 = z^N4
      cf = cf * z4
    end
  elseif n4 == 2
    if N % 8 != 0
      return
    end
    N8 = div(N, 8)
    z8 = z^N8
    cf = cf * (z8 - z8^3)
    nn = div(nn, 2)
    if mod(sign(n)*nn, 4) == 3
      cf = cf * (-z8^2)
    end
  elseif n4 == 3
    if n > 0
      if N % 4 != 0
        return
      end
      N4 = div(N, 4)
      z4 = z^N4
      cf = cf * (-z4)
    end
  end

  # Compute the coefficients of the Atlas irrationality 2*b_nn+1,
  # w.r.t. the N-th cyclotomic field.
  # (The underlying formula is due to a theorem of Gauss.)
  cfs = zeros(ZZRingElem, N)
  cfs[1] = 1
  q = div(N, nn)
  for k in 1:div(nn,2)
    pos = q * mod(k^2, nn) + 1
    cfs[pos] = cfs[pos] + 2
  end

  # Create the corresponding number field element.
  FF = parent(z.data)
  pol = FF.pol
  R = parent(pol)
  elm = mod(R(cfs), pol)

  return cf * QQAbElem(FF(elm), N)
end

"""
    quadratic_irrationality_info(a::QQAbModule.QQAbElem)

Return `(x, y, n)`, where `x`, `y` are of type `QQFieldElem` and `n` is
a squarefree integer, such that `a == x + y sqrt(n)` holds.

(We assume that the underlying primitive `N`-th root of unity that
is used to define `a` is identified with the complex number `exp(2*Pi*i/N)`,
where `i` is the imaginary unit.)
"""
function quadratic_irrationality_info(a::QQAbElem)
    n = a.c

    # Compute the Galois group generators of the `n`-th cyclotomic field.
    galgens = generators_galois_group_cyclotomic_field(n)

    # Start computing the orbit of `a` under the Galois group:
    # If the orbit length is larger than 2 then `a` is not in a
    # quadratic field, and we return nothing.
    cand = nothing
    for sigma in galgens
      img = a^sigma
      if img != a
        if isnothing(cand)
          cand = img
        elseif cand != img
          return
        end
      end
    end

    if cand === nothing
      # The value is rational.
      return (coeff(a.data, 0), 0, 1)
    end

    for sigma in galgens
      img = cand^sigma
      if img != a && img != cand
        # There are more than two Galois conjugates.
        return
      end
    end

    # We have a = x + y \sqrt{m} and cand = x - y \sqrt{m}.
    x = coeff(a.data + cand.data, 0) // 2
    root_multiple = a.data - x
    ysquarem = coeff(root_multiple^2, 0)  # QQFieldElem
    num = numerator(ysquarem)
    den = denominator(ysquarem)
    den_y = sqrt(den)
    m = sign(num)
    for (p, e) in collect(factor(num))
      if e % 2 == 1
        m = m * p
      end
    end
    y = sqrt(ysquarem // m)

    # It remains to compute the sign of y.
    # (This relies on the choice of the primitive n-th root of unity
    # as the complex number exp(2*pi*i/n).)
    y_std_m = y * square_root_in_cyclotomic_field(parent(a), Int(m), n).data
    if y_std_m != root_multiple
      @assert y_std_m == - root_multiple
      y = -y
    end

    return (x, y, m)
end

@doc raw"""
    reduce(val::QQAbElem, F::FinField)

Return the element of `F` that is the $p$-modular reduction of `val`,
where $p$ is the characteristic of `F`.
An exception is thrown if `val` cannot be reduced modulo $p$
or if the reduction does not lie in `F`.

# Examples
```jldocstring
julia> K, z = abelian_closure(QQ);

julia> F = GF(2, 3);

julia> reduce(z(7), F)
o
```
"""
function reduce(val::QQAbElem, F::FinField)
  p = characteristic(F)
  iso_0 = Oscar.iso_oscar_gap(parent(val))
  iso_p = Oscar.iso_oscar_gap(F)
  val_p = GAP.Globals.FrobeniusCharacterValue(iso_0(val), GAP.Obj(p))::GAP.Obj
  return preimage(iso_p, val_p)
end

#TODO: add reduction to alg. closure as soon as this is available!

end # module AbelianClosure

import .AbelianClosure:
       abelian_closure,
       QQAbAutomorphism,
       QQAbField,
       QQAbElem,
       set_variable!,
       get_variable

export abelian_closure
export get_variable
export QQAbAutomorphism
export QQAbElem
export QQAbField
export set_variable!
