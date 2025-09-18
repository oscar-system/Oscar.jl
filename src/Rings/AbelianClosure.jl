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
# saturations with allroots or allrootsNew (change this in the code)

abstract type CyclotomicField end

export CyclotomicField

module AbelianClosure 

using ..Oscar

import Base: +, *, -, //, ==, zero, one, ^, div, isone, iszero,
             deepcopy_internal, hash, reduce, isinteger

#import ..Oscar.AbstractAlgebra: promote_rule

import ..Oscar: AbstractAlgebra, add!, base_ring, base_ring_type, characteristic, elem_type, divexact, gen,
                has_preimage_with_preimage, is_root_of_unity, is_unit, mul!, neg!, parent,
                parent_type, promote_rule, root, root_of_unity, roots, @req

import Oscar: pretty, Lowercase

using Hecke
import Hecke: conductor, data, is_rational, is_integral, is_algebraic_integer

################################################################################
#
#  Types
#
################################################################################

@doc raw"""
    QQAbField

The type of the abelian closure of the rationals. An object of this type can
be constructed using [`abelian_closure(::QQField)`](@ref).
"""
@attributes mutable struct QQAbField{T} <: Nemo.Field # union of cyclotomic fields
  fields::Dict{Int, T} # Cache for the cyclotomic fields
  s::String

  function QQAbField{T}(fields::Dict{Int, T}) where T
    return new(fields)
  end
end

const _QQAb = QQAbField{AbsSimpleNumField}(Dict{Int, AbsSimpleNumField}())
const _QQAb_sparse = QQAbField{AbsNonSimpleNumField}(Dict{Int, AbsNonSimpleNumField}())

@doc raw"""
    QQAbFieldElem

Element type for the abelian closure of the rationals.
For more details see [`abelian_closure(::QQField)`](@ref).
"""
struct QQAbFieldElem{T} <: Nemo.FieldElem
  data::T                             # Element in cyclotomic field
  c::Int                              # Conductor of field
end
#T test that data really belongs to a cyclotomic field!

# This is a functor like object G with G(n) = primitive n-th root of unity

struct QQAbFieldGen{T}
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
    abelian_closure(QQ::QQField; sparse::Bool = false)

Return a pair `(K, z)` consisting of the abelian closure `K` of the rationals
and a generator `z` that can be used to construct primitive roots of unity in
`K`.
For each positive  integer `n`, `z(n)` is the primitive `n`-th root of unity
that corresponds to the complex number $\exp(2\pi i/n)$.
In particular, we have `z(n)^m = z(div(n, m))`, for all divisors `m` of `n`.

An optional keyword argument `sparse` can be set to `true` to switch to a 
sparse representation. Depending on the application this can be much faster
or slower.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
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

elem_type(::Type{QQAbField{AbsSimpleNumField}}) = QQAbFieldElem{AbsSimpleNumFieldElem}
parent_type(::Type{QQAbFieldElem{AbsSimpleNumFieldElem}}) = QQAbField{AbsSimpleNumField}
parent(::QQAbFieldElem{AbsSimpleNumFieldElem}) = _QQAb

elem_type(::Type{QQAbField{AbsNonSimpleNumField}}) = QQAbFieldElem{AbsNonSimpleNumFieldElem}
parent_type(::Type{QQAbFieldElem{AbsNonSimpleNumFieldElem}}) = QQAbField{AbsNonSimpleNumField}
parent(::QQAbFieldElem{AbsNonSimpleNumFieldElem}) = _QQAb_sparse

base_ring(::QQAbField) = Union{}
base_ring_type(::Type{<:QQAbField}) = typeof(Union{})

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

_variable(b::QQAbFieldElem{AbsSimpleNumFieldElem}) = Expr(:call, Symbol(_variable(_QQAb)), b.c)

function _variable(b::QQAbFieldElem{AbsNonSimpleNumFieldElem}) 
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

Hecke.data(a::QQAbFieldElem) = a.data

################################################################################
#
#  Creation of elements
#
################################################################################

@doc raw"""
    root_of_unity(K::QQAbField, n::Int)

Return the root of unity $\exp(2\pi i/n)$ as an element of `K`.

# Examples
```jldocstring
julia> K, z = abelian_closure(QQ);

julia> root_of_unity(K, 6) == z(6)
true
```
"""
function root_of_unity(K::QQAbField, n::Int)
  # Represent the root in the smallest possible cyclotomic field.
  if n % 2 == 0 && n % 4 != 0
    c = div(n, 2)
  else
    c = n
  end
  K, z = cyclotomic_field(K, c)
  if c != n
    z = -z^div(c+1, 2)
  end
  return QQAbFieldElem(z, c)
end

function root_of_unity2(K::QQAbField, c::Int)
  # This function returns the primitive root of unity e^(2*pi*i/c),
  # as an element of the `c`-th cyclotomic field,
  # also if `mod(c, 4 ) = 2` holds.
  K, z = cyclotomic_field(K, c)
  return QQAbFieldElem(z, c)
end

(z::QQAbFieldGen)(n::Int) = root_of_unity(z.K, n)
(z::QQAbFieldGen)(n::Int, r::Int) = z(n)^r

one(K::QQAbField) = K(1)

one(a::QQAbFieldElem) = one(parent(a))

function isone(a::QQAbFieldElem)
  return isone(data(a))
end

function iszero(a::QQAbFieldElem)
  return iszero(data(a))
end

zero(K::QQAbField) = K(0)

zero(a::QQAbFieldElem) = zero(parent(a))

function (K::QQAbField)(a::Union{ZZRingElem, QQFieldElem, Integer, Rational})
  return a*root_of_unity(K, 1)
end

function (K::QQAbField)(a::QQAbFieldElem)
  return a
end

(K::QQAbField)() = zero(K)

function (K::QQAbField)(a::AbsSimpleNumFieldElem)
  F = parent(a)

  # Cyclotomic fields are naturally embedded into `K`.
  fl, f = Hecke.is_cyclotomic_type(F)
  fl && return QQAbFieldElem(a, f)

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
    iszero(y) && return QQAbFieldElem(parent(a)(x), 1)
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
  print(pretty(io), "Sparse abelian closure of ", Lowercase(), QQ)
end

function Base.show(io::IO, a::QQAbField{AbsSimpleNumField})
  print(pretty(io), "Abelian closure of ", Lowercase(), QQ)
end

function Base.show(io::IO, a::QQAbFieldGen)
  print(pretty(io), "Generator of ", Lowercase(), a.K)
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

function AbstractAlgebra.expressify(b::QQAbFieldElem{AbsSimpleNumFieldElem}; context = nothing)
  a = data(b)
  return AbstractAlgebra.expressify(parent(parent(a).pol)(a), _variable(b), context = context)
end

function AbstractAlgebra.expressify(b::QQAbFieldElem{AbsNonSimpleNumFieldElem}; context = nothing)
  a = data(b)
  return AbstractAlgebra.expressify(a.data, _variable(b), context = context)
end

Oscar.@enable_all_show_via_expressify QQAbFieldElem

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

function coerce_up(K::AbsSimpleNumField, n::Int, a::QQAbFieldElem{AbsSimpleNumFieldElem})
  d = div(n, a.c)
  @assert n % a.c == 0
  #z_n^(d) = z_a
  R = parent(parent(data(a)).pol)
  return QQAbFieldElem{AbsSimpleNumFieldElem}(evaluate(R(data(a)), gen(K)^d), n)
end

function coerce_up(K::AbsNonSimpleNumField, n::Int, a::QQAbFieldElem{AbsNonSimpleNumFieldElem})
  d = div(n, a.c)
  @assert n % a.c == 0
  lk = get_attribute(parent(a.data), :decom)
  #gen(k, i) = gen(K, j)^n for the unique j s.th. gcd(lk[i], lK[j])
  # and n = lK[j]/lk[i]
  #z_n^(d) = z_a
  return QQAbFieldElem{AbsNonSimpleNumFieldElem}(evaluate(data(a).data, [ns_gen(K)^divexact(n, i) for i=lk]), n)
end


function coerce_down(K::AbsSimpleNumField, n::Int, a::QQAbFieldElem)
  throw(Hecke.NotImplemented())
end

function make_compatible(a::QQAbFieldElem{T}, b::QQAbFieldElem{T}) where {T}
  if a.c == b.c
    return a,b
  end
  d = lcm(a.c, b.c)
  K, = cyclotomic_field(parent(a), d)
  return coerce_up(K, d, a), coerce_up(K, d, b)
end

function minimize(::typeof(CyclotomicField), a::AbstractArray{AbsSimpleNumFieldElem})
  fl, c = Hecke.is_cyclotomic_type(parent(a[1]))
  @assert allequal(parent, a)
  @assert fl
  for p = prime_divisors(c)
    while c % p == 0
      K, _ = cyclotomic_field(Int(div(c, p)), cached = false)
      b = similar(a)
      OK = true
      for x = eachindex(a)
        y = Hecke.force_coerce_cyclo(K, a[x], Val(false))
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
  if is_conductor(c)
    return c
  else
    return div(c, 2)
  end
end

conductor(a::QQAbFieldElem) = conductor(data(a))

# What we want is the conductor of the domain of the map, but we need the map.
function conductor(phi::MapFromFunc{T, QQAbField{T}}) where T
  return lcm([conductor(phi(x)) for x in gens(domain(phi))])
end

################################################################################
#
#  Conversions to `ZZRingElem` and `QQFieldElem` (like for `AbsSimpleNumFieldElem`)
#
################################################################################

(R::QQField)(a::QQAbFieldElem) = R(a.data)
(R::ZZRing)(a::QQAbFieldElem) = R(a.data)

################################################################################
#
#  Conversion to `QQBarFieldElem`
#
#  We assume the natural embedding of cyclotomic fields into `QQBarField()`
#  that is given by the (documented) fact that
#  `root_of_unity(F::QQBarField, n::Int)` is $\exp(2 \pi i / n)$.
#
################################################################################

function (F::QQBarField)(a::QQAbFieldElem)
  N = a.c
  cfs = Oscar.coefficients(a.data)
  r = root_of_unity(F, N)
  pow = one(F)
  tmp = zero(F)
  res = cfs[1] * pow
  for i in 2:length(cfs)
    pow = mul!(pow, r)
    res = addmul!(res, cfs[i], pow, tmp)
  end
  return res
end

Nemo.QQBarFieldElem(a::QQAbFieldElem) = algebraic_closure(QQ)(a)

################################################################################
#
#  Conversion to `Float64`, `ComplexF64`
#
#  We first convert to a `QQBarFieldElem` and then use that it supports the
#  conversions in question.
#

Core.Float64(a::QQAbFieldElem) = Float64(QQBarFieldElem(a))
Base.ComplexF64(a::QQAbFieldElem) = ComplexF64(QQBarFieldElem(a))

################################################################################
#
#  Ring interface functions
#
################################################################################

is_unit(a::QQAbFieldElem) = !iszero(a)

canonical_unit(a::QQAbFieldElem) = a

################################################################################
#
#  Minimal polynomial
#
################################################################################

Hecke.minpoly(a::QQAbFieldElem) = minpoly(data(a))

################################################################################
#
#  Subfields and such
#
################################################################################
"""
To parametrize subfields of the n-cyclotomic field we use blocks:
 - pick a prime n = 1 % p (so thre is a n-th root of 1 mod p)
 - pick an n-th root Z of one
 - the roots are exactly Z^i for i coprime to n as the images
   of the n-th root of 1 in the QQab field
This way the ordering of the roots is mentained across different primes
(choice of Z correpsonds to fixing a prime ideal above p)
"""
struct RootData
   p::Int # n = 1 % p
   r::Vector{Int} # the roots in this order as described above

   function RootData(n::Int, p::Int)
     @assert p % n == 1
     k = Native.GF(p)
     @assert((p-1)%n == 0)
     lf = prime_divisors(n)
     local Z, Zn
     while true
       Z = rand(k)
       iszero(Z) && continue
       Zn = Z^divexact(p-1, n)
       if all(x->!isone(Zn^divexact(n, x)), lf)
         break
       end
     end
     r, mr = quo(ZZ, n)
     u, mu = unit_group(r)
     c = sort(Int[preimage(mr, mu(g)) for g = u])
     
     z = Int(lift(Zn))
     return new(p, Int[powermod(z, x, p) for x = c])
   end
end

function block_system(a::AbsSimpleNumFieldElem, rd::RootData)
  p = rd.p
  if denominator(a) % p == 0
    return Vector{Vector{Int}}()
  end
  kp = Native.GF(p; check = false, cached = false)
  Qx = parent(defining_polynomial(parent(a)))
  pol = Qx(a)
  v = [map_coefficients(kp, pol)(t) for t = rd.r]
  D = Dict{fpFieldElem, Vector{Int}}()
  for i=1:length(v)
    if haskey(D, v[i])
      push!(D[v[i]], i)
    else
      D[v[i]] = [i]
    end
  end
  s = sort(collect(values(D)), lt = (a,b) -> isless(a[1], b[1]))
  if any(x->length(x) != length(s[1]), s)
    return Vector{Vector{Int}}()
  end
  return s
end

function block_system(k::AbsSimpleNumField, a::QQAbFieldElem)
  if degree(k) == 1
    return [[1]]
  end
  A = parent(a)
  rda = get_attribute(A, :RootData)
  if rda === nothing
    rda = Dict{Int, Vector{RootData}}()
    set_attribute!(A, :RootData => rda)
  end
  n = Hecke.is_cyclotomic_type(parent(k(data(a))))[2]
  if !haskey(rda, n)
    rda[n] = Vector{Vector{Int}}()
  end
  p = 1
  for rd = rda[n]
    p = max(p, rd.p)
    b = block_system(k(data(a)), rd)
    if length(b) != 0
      return b
    end
  end
  nq = 1
  for q = PrimesSet(p+1, -1, n, 1)
    rd = RootData(n, q)
    push!(rda[n], rd)
    b = block_system(k(data(a)), rd)
    if length(b) != 0
      return b
    end
    nq += 1
    if nq > 100 # not plausible, s.th. is going badly wrong
      #a will need to be divisible at 100 primes...
      error("dnw")
    end
  end
end

function intersect_block_systems(a::Vector{Vector{Int}}, b::Vector{Vector{Int}})
  c = [intersect(x, y) for x = a for y = b]
  return sort([x for x = c if length(x) > 0], lt = (a,b) -> isless(a[1], b[1]))
end

function Oscar.sub(K::QQAbField, s::Vector{<:QQAbFieldElem}; cached::Bool = true)
  f = lcm([Hecke.is_cyclotomic_type(parent(data(x)))[2] for x = s])
  k = cyclotomic_field(f)[1]
  b = [[i for i = 1:degree(k)]] #block system for QQ as a subfield
  pe = zero(K)
  for mu = s
    bs = block_system(k, mu)
    if issubset(b[1], bs[1])
      continue #nothing new in element
    end
    b = intersect_block_systems(b, bs)
    i = 1
    while block_system(k, pe+i*mu) != b
      i += 1
      if i> 10 #in theory the number of failures is finite
          #this is just a safety valve, could be removed
        error("dnw")
      end
    end
    pe += i*mu
  end
  if iszero(pe) #to catch QQ as a subfield, we prefer x-1 over x
    pe = one(K)
  end
  if cached
    old = get_attribute(K, :subfields)
    if old === nothing
      old = Dict{Tuple{Int, Vector{Int}}, Map}()
      set_attribute!(K, :subfields=>old)
    end
    if haskey(old, (f, b[1]))
      hh = old[(f, b[1])]
      return domain(hh), hh
    end
  end
  g = minpoly(pe)
  @assert degree(g) == length(b)
  s, _ = number_field(g; check = false, cached = false)
  h = hom(s, k, k(pe.data))
  hh = MapFromFunc(s, K, x->K(h(x)), y-> preimage(h, k(y)))
  if cached
    old = get_attribute(K, :subfields)
    old[(f, b[1])] = hh
  end
  return s, hh
end
    
################################################################################
#
#  Syntactic sugar
#
################################################################################

function Hecke.number_field(K::QQField, a::QQAbFieldElem; cached::Bool = false)
  return number_field(K, [a]; cached)
end

function Hecke.number_field(K::QQField, a::AbstractVector{<: QQAbFieldElem}; cached::Bool = false)
  k, mp = sub(parent(a[1]), a; cached)
  return k, gen(k)
end

Base.getindex(::QQField, a::QQAbFieldElem) = number_field(QQ, a)
Base.getindex(::QQField, a::Vector{QQAbFieldElem{T}}) where T = number_field(QQ, a)
Base.getindex(::QQField, a::QQAbFieldElem...) = number_field(QQ, [x for x in a])

################################################################################
#
#   Unary operators
#
################################################################################

function -(a::QQAbFieldElem)
  return QQAbFieldElem(-data(a), a.c)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(data(a) + data(b), a.c)
end

function -(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(data(a) - data(b), a.c)
end

function *(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(data(a) * data(b), a.c)
end

function ^(a::QQAbFieldElem, n::Integer)
  return QQAbFieldElem(data(a)^n, a.c)
end

function ^(a::QQAbFieldElem, n::ZZRingElem)
  return a^Int(n)
end

################################################################################
#
#   Exact division
#
################################################################################

function //(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(a.data//b.data, a.c)
end

function div(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(a.data//b.data, a.c)
end

function divexact(a::QQAbFieldElem, b::QQAbFieldElem; check::Bool = true)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(divexact(a.data, b.data), a.c)
end

function inv(a::QQAbFieldElem)
  return QQAbFieldElem(inv(data(a)), a.c)
end

################################################################################
#
#  Unsafe operations
#
################################################################################

function neg!(a::QQAbFieldElem)
  return QQAbFieldElem(neg!(a.data), a.c)
end

function add!(c::QQAbFieldElem, a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  if c.c != a.c
    return a + b
  else
    return QQAbFieldElem(add!(c.data, a.data, b.data), a.c)
  end
end

function add!(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(add!(a.data, b.data), a.c)
end

function sub!(c::QQAbFieldElem, a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  if c.c != a.c
    return a - b
  else
    return QQAbFieldElem(sub!(c.data, a.data, b.data), a.c)
  end
end

function sub!(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(sub!(a.data, b.data), a.c)
end

function mul!(c::QQAbFieldElem, a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  if c.c != a.c
    return a * b
  else
    return QQAbFieldElem(mul!(c.data, a.data, b.data), a.c)
  end
end

function mul!(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return QQAbFieldElem(mul!(a.data, b.data), a.c)
end

################################################################################
#
#  Ad hoc binary operations
#
################################################################################

for T in (ZZRingElem, QQFieldElem, Integer, Rational)
  @eval begin
    +(a::QQAbFieldElem, b::$T) = QQAbFieldElem(data(a)+b, a.c)
    +(a::$T, b::QQAbFieldElem) = b+a

    -(a::QQAbFieldElem, b::$T) = QQAbFieldElem(data(a)-b, a.c)
    -(a::$T, b::QQAbFieldElem) = QQAbFieldElem(a-data(b), b.c)

    *(a::QQAbFieldElem, b::$T) = QQAbFieldElem(data(a)*b, a.c)
    *(a::$T, b::QQAbFieldElem) = b*a

    //(a::QQAbFieldElem, b::$T) = QQAbFieldElem(data(a)//b, a.c)
    //(a::$T, b::QQAbFieldElem) = QQAbFieldElem(a//data(b), b.c)

    divexact(a::QQAbFieldElem, b::$T; check::Bool = true) = QQAbFieldElem(data(a)/b, a.c)
    divexact(a::$T, b::QQAbFieldElem; check::Bool = true) = QQAbFieldElem(a/data(b), b.c)
  end
end

################################################################################
#
#  Comparison
#
################################################################################

function ==(a::QQAbFieldElem, b::QQAbFieldElem)
  a, b = make_compatible(a, b)
  return a.data == b.data
end

function ==(a::QQAbFieldElem, b::Union{ZZRingElem, QQFieldElem, Integer, Rational})
  return data(a) == b
end

function ==(a::Union{ZZRingElem, QQFieldElem, Integer, Rational}, b::QQAbFieldElem)
  return b == a
end

hash(a::QQAbFieldElem, h::UInt) = hash(minimize(CyclotomicField, a.data), h)

################################################################################
#
#  Copy and deepcopy
#
################################################################################

function Base.copy(a::QQAbFieldElem)
  return QQAbFieldElem(data(a), a.c)
end

function Base.deepcopy_internal(a::QQAbFieldElem, dict::IdDict)
  return QQAbFieldElem(deepcopy_internal(data(a), dict), a.c)
end

################################################################################
#
#  Promotion rules
#
################################################################################

AbstractAlgebra.promote_rule(::Type{QQAbFieldElem}, ::Type{Int}) = QQAbFieldElem

AbstractAlgebra.promote_rule(::Type{QQAbFieldElem}, ::Type{ZZRingElem}) = QQAbFieldElem

AbstractAlgebra.promote_rule(::Type{QQAbFieldElem}, ::Type{QQFieldElem}) = QQAbFieldElem

###############################################################################
#
#  isinteger, is_rational, is_integral
#
###############################################################################

@doc raw"""
    isinteger(a::QQAbFieldElem)

Return whether $a$ is an integer.
"""
isinteger(a::QQAbFieldElem) = isinteger(data(a))

@doc raw"""
    is_rational(a::QQAbFieldElem)

Return whether $a$ is a rational number.
"""
is_rational(a::QQAbFieldElem) = is_rational(data(a))

@doc raw"""
    is_integral(a::QQAbFieldElem)

Returns whether $a$ is integral, that is, whether the minimal
polynomial of $a$ has integral coefficients.
"""
is_integral(a::QQAbFieldElem) = is_integral(data(a))

@doc raw"""
    is_algebraic_integer(a::QQAbFieldElem)

Return whether $a$ is an algebraic integer.
"""
is_algebraic_integer(a::QQAbFieldElem) = is_integral(a)

###############################################################################
#
#   Functions for computing roots
#
###############################################################################

function Oscar.root(a::QQAbFieldElem, n::Int)
  @req is_root_of_unity(a) "Element must be a root of unity"
  o = Oscar.order(a)
  l = o*n
  mu = root_of_unity2(parent(a), Int(l))
  return mu
end

function Oscar.roots(f::PolyRingElem{QQAbFieldElem{T}}) where T
  QQAb = base_ring(f)
  c = reduce(lcm, map(conductor, AbstractAlgebra.coefficients(f)), init = Int(1))
  k, z = cyclotomic_field(QQAb, c)
  f = map_coefficients(x->k(x.data), f)
  lf = factor(f)
  #we need to find the correct cyclotomic field...
  #can't use ray_class_group in k as this is expensive (needs class group)
  #need absolute norm group
  QQ = rationals_as_number_field()[1]

  C = cyclotomic_field(ClassField, c)

  rts = QQAbFieldElem{T}[]

  for (g, _) in lf
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
        q, mqq = quo(q, [degree(x)*mq(P) for (x, _) in l], false)
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
  return rts::Vector{QQAbFieldElem{T}}
end

function Oscar.roots(a::QQAbFieldElem{T}, n::Int) where {T}
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
    _, x = polynomial_ring(parent(a); cached = false)
    fl || return roots(x^n-a)::Vector{QQAbFieldElem{T}}
    b = gens(Hecke.inv(i))[end]
    c = QQAbFieldElem(b, a.c)
    corr = Hecke.inv(c)
    a *= c^n
    fl = is_root_of_unity(a)
    fl || return (corr .* roots(x^n-a))::Vector{QQAbFieldElem{T}}
  end
  
  o = order(a)
  l = o*n
  mu = root_of_unity(parent(a), Int(l))
  A = QQAbFieldElem[]
  if l==1 && mu==a
    push!(A, mu)
  end
  for k = 0:(l-1)
    el = mu^k
    if el^n == a
      push!(A, el)
    end
  end
  return [x*corr for x = A]::Vector{QQAbFieldElem{T}}
end

function is_root_of_unity(a::QQAbFieldElem)
  return is_torsion_unit(a.data, true)
  #=
  b = a^a.c
  return b.data == 1 || b.data == -1
  =#
end

function Oscar.order(a::QQAbFieldElem)
  f = Nemo.factor(ZZRingElem(2*a.c))
  o = 1
  for (p, e) in f
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
# which is already implemented for QQAbFieldElem directly

function Oscar.sqrt(a::QQAbFieldElem)
  sqrt = Oscar.roots(a, 2)
  if is_empty(sqrt)
    error("Element $a does not have a square root")
  end
  return sqrt[1]
end

function Oscar.cbrt(a::QQAbFieldElem)
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
# If `F` has conductor `N` then assume that `x.c == N` holds.
# If `F` is a cyclotomic field with conductor `N` then assume that
# `x == QQAbFieldElem(gen(F), N)`.
# (Use that the powers of this element form a basis of the field.)
function _embedding(F::QQField, K::QQAbField{AbsSimpleNumField},
                    x::QQAbFieldElem{AbsSimpleNumFieldElem})
  C1, _ = cyclotomic_field(1)

  f = function(x::QQFieldElem)
    return QQAbFieldElem(C1(x), 1)
  end

  finv = function(x::QQAbFieldElem; throw_error::Bool = true)
    res = Hecke.force_coerce_cyclo(C1, data(x), Val(false))
    throw_error && res === nothing && error("element has no preimage")
    return res
  end

  return MapFromFunc(F, K, f, finv)
end

function _embedding(F::AbsSimpleNumField, K::QQAbField{AbsSimpleNumField},
                    x::QQAbFieldElem{AbsSimpleNumFieldElem})
  fl, n = Hecke.is_cyclotomic_type(F)
  if fl
    # This is cheaper.
    f = function(x::AbsSimpleNumFieldElem)
      return QQAbFieldElem(x, n)
    end

    finv = function(x::QQAbFieldElem; throw_error::Bool = true)
      res = Hecke.force_coerce_cyclo(F, data(x), Val(false))
      throw_error && res === nothing && error("element has no preimage")
      return res
    end
  else
    # `F` is expected to be a proper subfield of a cyclotomic field.
    n = x.c
    x = data(x)
    Kn, = AbelianClosure.cyclotomic_field(K, n)
    powers = [Hecke.coefficients(Hecke.force_coerce_cyclo(Kn, x^i))
              for i in 0:degree(F)-1]
    c = transpose(matrix(QQ, powers))
    R = parent(F.pol)

    f = function(z::AbsSimpleNumFieldElem)
      return QQAbFieldElem(evaluate(R(z), x), n)
    end

    finv = function(x::QQAbFieldElem; throw_error::Bool = true)
      # Write `x` w.r.t. the n-th cyclotomic field ...
      g = gcd(x.c, n)
      Kg, = AbelianClosure.cyclotomic_field(K, g)
      x = Hecke.force_coerce_cyclo(Kg, data(x), Val(false))
      if x === nothing
        throw_error && error("element has no preimage")
        return
      end
      x = Hecke.force_coerce_cyclo(Kn, x)
      # ... and then w.r.t. `F`
      a = Hecke.coefficients(x)
      fl, sol = can_solve_with_solution(c, matrix(QQ, length(a), 1, a); side = :right)
      if !fl
        throw_error && error("element has no preimage")
        return
      end
      b = transpose(sol)
      b = [b[i] for i in 1:length(b)]
      return F(b)
    end
  end
  return MapFromFunc(F, K, f, finv)
end

# The following works only if `mp.g` admits a second argument,
# which is the case if `mp` has been constructed by `_embedding` above.
function has_preimage_with_preimage(mp::MapFromFunc{AbsSimpleNumField, QQAbField{AbsSimpleNumField}}, x::QQAbFieldElem{AbsSimpleNumFieldElem})
  pre = mp.g(x, throw_error = false)
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

struct QQAbAutomorphism
  exp::Int
end

Oscar.hom(K::QQAbField, L::QQAbField, k::Int) = QQAbAutomorphism(k)

(f::QQAbAutomorphism)(a::QQAbFieldElem) = a^f

function ^(val::QQAbFieldElem, sigma::QQAbAutomorphism)
  k = sigma.exp
  n = val.c
  g = gcd(k, n)
  if g != 1
    # Replace `k` by an equivalent one that is coprime to `n`.
    n0 = 1
    n1 = n
    for (p, exp) in Oscar.factor(g)
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
  return QQAbFieldElem(F(R(res)), n)
end

Base.conj(elm::QQAbFieldElem) = elm^QQAbAutomorphism(-1)

Base.isreal(elm::QQAbFieldElem) = conj(elm) == elm

# compare real `QQAbFieldElem`s
function Base.isless(a::QQAbFieldElem, b::QQAbFieldElem)
  F = QQBarField()
  return Base.isless(F(a), F(b))
end

_isless_via_qqbar(a, b) = Base.isless(QQBarFieldElem(a), QQBarFieldElem(b))

for T in (QQFieldElem, ZZRingElem, Int, Integer, Rational)
  @eval begin
    Base.isless(a::QQAbFieldElem, b::$T) = _isless_via_qqbar(a, b)
    Base.isless(a::$T, b::QQAbFieldElem) = _isless_via_qqbar(a, b)
  end
end

AbstractAlgebra.is_positive(a::QQAbFieldElem) = _isless_via_qqbar(0, a)
AbstractAlgebra.is_negative(a::QQAbFieldElem) = _isless_via_qqbar(a, 0)

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
  for (p,e) in factor(n)
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
  # (The underlying formula is due to a theorem of Gauss,
  # see Chapter IV, § 3, QS 4 in [Lan70](@cite).)
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

  return cf * QQAbFieldElem(FF(elm), N)
end

"""
    quadratic_irrationality_info(a::QQAbFieldElem)

Return `(x, y, n)`, where `x`, `y` are of type `QQFieldElem` and `n` is
a squarefree integer, such that `a == x + y sqrt(n)` holds.

(We assume that the underlying primitive `N`-th root of unity that
is used to define `a` is identified with the complex number `exp(2*Pi*i/N)`,
where `i` is the imaginary unit.)
"""
function quadratic_irrationality_info(a::QQAbFieldElem)
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
    for (p, e) in factor(num)
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
    reduce(val::QQAbFieldElem, F::FinField)

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
function reduce(val::QQAbFieldElem, F::FinField)
  p = characteristic(F)
  iso_0 = Oscar.iso_oscar_gap(parent(val))
  iso_p = Oscar.iso_oscar_gap(F)
  val_p = GAP.Globals.FrobeniusCharacterValue(iso_0(val), GAP.Obj(p))::GAP.Obj
  return preimage(iso_p, val_p)
end

#TODO: add reduction to alg. closure as soon as this is available!

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(K::QQAbField)
  ns = rand(1:8, 3)
  zs = map(n -> sum(rand(-10:10) * gen(K)(n)^rand(1:n) for j in 1:10), ns)
  return sum(zs)
end

end # module AbelianClosure

import .AbelianClosure:
       abelian_closure,
       QQAbAutomorphism,
       QQAbField,
       QQAbFieldElem,
       set_variable!,
       get_variable

export abelian_closure
export get_variable
export QQAbAutomorphism
export QQAbFieldElem
export QQAbField
export set_variable!
