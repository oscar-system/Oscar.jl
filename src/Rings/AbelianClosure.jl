###############################################################################
#
#  Abelian closure of the rationals
#
###############################################################################

# This is an implementation of Q^ab, the abelian closure of the rationals,
# which is modelled as the union of cyclotomic fields.
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

module AbelianClosure 

using ..Oscar, Markdown

import Base: +, *, -, //, ==, zero, one, ^, div, isone, iszero, deepcopy_internal, hash

#import ..Oscar.AbstractAlgebra: promote_rule

import ..Oscar: addeq!, is_unit, parent_type, elem_type, gen, root_of_unity,
                root, divexact, mul!, roots, is_root_of_unity, promote_rule,
                AbstractAlgebra, parent
using Hecke
import Hecke: conductor, data

################################################################################
#
#  Types
#
################################################################################

@attributes mutable struct QQAbField{T} <: Nemo.Field # union of cyclotomic fields
  s::String
  fields::Dict{Int, T} # Cache for the cyclotomic fields

  function QQAbField{T}(s::String, fields::Dict{Int, T}) where T
    return new(s, fields)
  end
end

const _QQAb = QQAbField{AnticNumberField}("ζ", Dict{Int, AnticNumberField}())
const _QQAb_sparse = QQAbField{NfAbsNS}("ζ", Dict{Int, NfAbsNS}())

mutable struct QQAbElem{T} <: Nemo.FieldElem
  data::T                             # Element in cyclotomic field
  c::Int                              # Conductor of field
end

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

@doc Markdown.doc"""
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
ζ(36)

julia> K, z = abelian_closure(QQ, sparse = true);

julia> z(36)
-ζ(36, 9)*ζ(36, 4)^4 - ζ(36, 9)*ζ(36, 4)

```
"""
function abelian_closure(::FlintRationalField; sparse::Bool = false) 
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
gen(K::QQAbField{AnticNumberField}) = _QQAbGen
gen(K::QQAbField{NfAbsNS}) = _QQAbGen_sparse

"""
    gen(K::QQAbField, s::String)

Return the generator of the abelian closure `K` that can be used to construct
primitive roots of unity. The string `s` will be used during printing.
"""
function gen(K::QQAbField, s::String)
  K.s = s
  return gen(K)
end

################################################################################
#
#  Parent and element functions
#
################################################################################

elem_type(::Type{QQAbField{AnticNumberField}}) = QQAbElem{nf_elem}
elem_type(::QQAbField{AnticNumberField}) = QQAbElem{nf_elem}
parent_type(::Type{QQAbElem{nf_elem}}) = QQAbField{AnticNumberField}
parent_type(::QQAbElem{nf_elem}) = QQAbField{AnticNumberField}
parent(::QQAbElem{nf_elem}) = _QQAb

elem_type(::Type{QQAbField{NfAbsNS}}) = QQAbElem{NfAbsNSElem}
elem_type(::QQAbField{NfAbsNS}) = QQAbElem{NfAbsNSElem}
parent_type(::Type{QQAbElem{NfAbsNSElem}}) = QQAbField{NfAbsNS}
parent_type(::QQAbElem{NfAbsNSElem}) = QQAbField{NfAbsNS}
parent(::QQAbElem{NfAbsNSElem}) = _QQAb_sparse

################################################################################
#
#  Field access
#
################################################################################

_variable(K::QQAbField) = K.s
_variable(b::QQAbElem{nf_elem}) = Expr(:call, Symbol(_variable(_QQAb)), b.c)

function _variable(b::QQAbElem{NfAbsNSElem}) 
  k = parent(b.data)
  lc = get_attribute(k, :decom)
  n = get_attribute(k, :cyclo)
  return [Expr(:call, Symbol(_variable(_QQAb)), n, divexact(n, i)) for i = lc]
end

function Hecke.cyclotomic_field(K::QQAbField{AnticNumberField}, c::Int)
  if haskey(K.fields, c)
    k = K.fields[c]
    return k, gen(k)
  else
    k, z = CyclotomicField(c, string(K.s, "(", c, ")"), cached = false)
    K.fields[c] = k
    return k, z
  end
end

function ns_gen(K::NfAbsNS)
  #z_pq^p = z_q and z_pg^q = z_p
  #thus z_pq = z_p^a z_q^b implies
  #z_pq^p = z_q^pb, so pb = 1 mod q
  #so:
  lc = get_attribute(K, :decom)
  n = get_attribute(K, :cyclo)
  return prod(gen(K, i)^invmod(divexact(n, lc[i]), lc[i]) for i=1:length(lc))
end

function Hecke.cyclotomic_field(K::QQAbField{NfAbsNS}, c::Int)
  if haskey(K.fields, c)
    k = K.fields[c]
    return k, ns_gen(k)
  else
    k, _ = cyclotomic_field(NonSimpleNumField, c, string(K.s))
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

function root_of_unity(K::QQAbField{AnticNumberField}, n::Int)
  if n % 2 == 0 && n % 4 != 0
    c = div(n, 2)
  else
    c = n
  end
  K, z = cyclotomic_field(K, c)
  if c == n
    return QQAbElem{nf_elem}(z, c)
  else
    return QQAbElem{nf_elem}(-z, c)
  end
end

function root_of_unity(K::QQAbField{NfAbsNS}, n::Int)
  if n % 2 == 0 && n % 4 != 0
    c = div(n, 2)
  else
    c = n
  end
  K, z = cyclotomic_field(K, c)
  if c == n
    return QQAbElem{NfAbsNSElem}(z, c)
  else
    return QQAbElem{NfAbsNSElem}(-z, c)
  end
end


function root_of_unity2(K::QQAbField, c::Int)
  # This function returns the primitive root of unity e^(2*pi*i/n)
  K, z = cyclotomic_field(K, c)
  return QQAbElem(z, c)
end

(z::QQAbFieldGen)(n::Int) = root_of_unity(z.K, n)
(z::QQAbFieldGen)(n::Int, r::Int) = z(n)^r

one(K::QQAbField) = _QQAb(1)

one(a::QQAbElem) = one(parent(a))

function isone(a::QQAbElem)
  return isone(data(a))
end

function iszero(a::QQAbElem)
  return iszero(data(a))
end

zero(K::QQAbField) = _QQAb(0)

zero(a::QQAbElem) = zero(parent(a))

function (K::QQAbField)(a::Union{fmpz, fmpq, Integer, Rational})
  return a*root_of_unity(K, 1)
end

function (K::QQAbField)(a::QQAbElem)
  return a
end

(K::QQAbField)() = zero(K)

################################################################################
#
#  String I/O
#
################################################################################

function Base.show(io::IO, a::QQAbField{NfAbsNS})
  print(io, "(Sparse) abelian closure of Q")
end
function Base.show(io::IO, a::QQAbField{AnticNumberField})
  print(io, "Abelian closure of Q")
end

function Base.show(io::IO, a::QQAbFieldGen)
  if isa(a.K, QQAbField{AnticNumberField})
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
  ss = K.s
  K.s = s
  return ss
end

"""
    get_variable(K::QQAbField)

Return the string used to print the primitive n-th root of the abelian closure
of the rationals.
"""
get_variable(K::QQAbField) = _variable(K)

function AbstractAlgebra.expressify(b::QQAbElem{nf_elem}; context = nothing)
  a = data(b)
  return AbstractAlgebra.expressify(parent(parent(a).pol)(a), _variable(b), context = context)
end

function AbstractAlgebra.expressify(b::QQAbElem{NfAbsNSElem}; context = nothing)
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

function coerce_up(K::AnticNumberField, n::Int, a::QQAbElem{nf_elem})
  d = div(n, a.c)
  @assert n % a.c == 0
  #z_n^(d) = z_a
  R = parent(parent(data(a)).pol)
  return QQAbElem{nf_elem}(evaluate(R(data(a)), gen(K)^d), n)
end

function coerce_up(K::NfAbsNS, n::Int, a::QQAbElem{NfAbsNSElem})
  d = div(n, a.c)
  @assert n % a.c == 0
  lk = get_attribute(parent(a.data), :decom)
  #gen(k, i) = gen(K, j)^n for the unique j s.th. gcd(lk[i], lK[j])
  # and n = lK[j]/lk[i]
  #z_n^(d) = z_a
  return QQAbElem{NfAbsNSElem}(evaluate(data(a).data, [ns_gen(K)^divexact(n, i) for i=lk]), n)
end


function coerce_down(K::AnticNumberField, n::Int, a::QQAbElem)
  throw(Hecke.NotImplemented())
end

function make_compatible(a::QQAbElem, b::QQAbElem)
  if a.c == b.c
    return a,b
  end
  d = lcm(a.c, b.c)
  K, = cyclotomic_field(parent(a), d)
  return coerce_up(K, d, a), coerce_up(K, d, b)
end

function minimize(::typeof(CyclotomicField), a::AbstractArray{nf_elem})
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

function minimize(::typeof(CyclotomicField), a::MatElem{nf_elem})
  return matrix(minimize(CyclotomicField, a.entries))
end

function minimize(::typeof(CyclotomicField), a::nf_elem)
  return minimize(CyclotomicField, [a])[1]
end

conductor(a::nf_elem) = conductor(parent(minimize(CyclotomicField, a)))

function conductor(k::AnticNumberField)
  f, c = Hecke.is_cyclotomic_type(k)
  f || error("field is not of cyclotomic type")
  return c
end

conductor(a::QQAbElem) = conductor(data(a))

################################################################################
#
#  Conversions to `fmpz` and `fmpq` (like for `nf_elem`)
#
################################################################################

(R::FlintRationalField)(a::QQAbElem) = R(a.data)
(R::FlintIntegerRing)(a::QQAbElem) = R(a.data)


################################################################################
#
#  Ring interface functions
#
################################################################################

is_unit(a::QQAbElem) = !iszero(a)

canonical_unit(a::QQAbElem) = a

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

function ^(a::QQAbElem, n::fmpz)
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

*(a::fmpz, b::QQAbElem) = QQAbElem(b.data*a, b.c)

*(a::fmpq, b::QQAbElem) = QQAbElem(b.data*a, b.c)

*(a::Integer, b::QQAbElem) = QQAbElem(data(b) * a, b.c)

*(a::Rational, b::QQAbElem) = QQAbElem(data(b) * a, b.c)

*(a::QQAbElem, b::fmpz) = b*a

*(a::QQAbElem, b::fmpq) = b*a

*(a::QQAbElem, b::Integer) = b*a

*(a::QQAbElem, b::Rational) = b*a

+(a::fmpz, b::QQAbElem) = QQAbElem(b.data + a, b.c)

+(a::fmpq, b::QQAbElem) = QQAbElem(b.data + a, b.c)

+(a::Integer, b::QQAbElem) = QQAbElem(data(b) + a, b.c)

+(a::Rational, b::QQAbElem) = QQAbElem(data(b) + a, b.c)

+(a::QQAbElem, b::fmpz) = b + a

+(a::QQAbElem, b::fmpq) = b + a

+(a::QQAbElem, b::Integer) = b + a

+(a::QQAbElem, b::Rational) = b + a

-(a::fmpz, b::QQAbElem) = QQAbElem(-(a, data(b)), b.c)

-(a::fmpq, b::QQAbElem) = QQAbElem(-(a, data(b)), b.c)

-(a::Integer, b::QQAbElem) = QQAbElem(-(a, data(b)), b.c)

-(a::Rational, b::QQAbElem) = QQAbElem(-(a, data(b)), b.c)

-(a::QQAbElem, b::fmpz) = QQAbElem(-(data(a), b), a.c)

-(a::QQAbElem, b::fmpq) = QQAbElem(-(data(a), b), a.c)

-(a::QQAbElem, b::Integer) = QQAbElem(-(data(a), b), a.c)

-(a::QQAbElem, b::Rational) = QQAbElem(-(data(a), b), a.c)

//(a::QQAbElem, b::fmpz) = QQAbElem(data(a)//b, a.c)

//(a::QQAbElem, b::fmpq) = QQAbElem(data(a)//b, a.c)

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

function ==(a::QQAbElem, b::Union{fmpz, fmpq, Integer, Rational})
  return data(a) == b
end

function ==(a::Union{fmpz, fmpq, Integer, Rational}, b::QQAbElem)
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

AbstractAlgebra.promote_rule(::Type{QQAbElem}, ::Type{fmpz}) = QQAbElem

AbstractAlgebra.promote_rule(::Type{QQAbElem}, ::Type{fmpq}) = QQAbElem

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

function Oscar.roots(f::PolyElem{QQAbElem{T}}) where T
  QQAb = base_ring(f)
  c = reduce(lcm, map(conductor, coefficients(f)), init = Int(1))
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
    c = reduce(lcm, map(conductor, coefficients(g)), init = Int(1))
    #so THIS factor lives in cyclo(c)
    k, z = cyclotomic_field(QQAb, c)
    d = numerator(norm(k(discriminant(g))))

    R, mR = ray_class_group(lcm(d, c)*maximal_order(QQ), infinite_places(QQ), 
                                        n_quo = degree(g)*degree(k))
    q, mq = quo(R, [R[0]])
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
        q, mqq = quo(q, [degree(x)*mq(P) for x = keys(l.fac)])
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
    _, x = PolynomialRing(parent(a), cached = false)
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
  f = Nemo.factor(fmpz(2*a.c))
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
  data = val.data  # nf_elem
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
  z = root_of_unity(F, N)
  cf = 1
  sqf = 1
  for (p,e) in collect(factor(n))
    cf = cf * p^div(e, 2)
    if e % 2 != 0
      sqf = sqf * p
    end
  end
  nn = Int(sqf)  # nn is positive and squarefree
  @assert N % nn == 0

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
    if n < 0
      cf = cf * z8^2
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
  cfs = zeros(fmpz, N)
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

Return `(x, y, n)`, where `x`, `y` are of type `fmpq` and `n` is
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

    if cand == nothing
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
    ysquarem = coeff(root_multiple^2, 0)  # fmpq
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

end # module AbelianClosure

import .AbelianClosure:
       abelian_closure,
       QQAbAutomorphism,
       QQAbField,
       QQAbElem,
       set_variable!,
       get_variable

export abelian_closure,
       QQAbAutomorphism,
       QQAbField,
       QQAbElem,
       set_variable!,
       get_variable
