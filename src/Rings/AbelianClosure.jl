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

using ..Oscar

import Base: +, *, -, //, ==, zero, one, ^, div, isone, iszero, deepcopy_internal

#import ..Oscar.AbstractAlgebra: promote_rule

import ..Oscar: addeq!, isunit, parent_type, elem_type, gen, root_of_unity,
                root, divexact, mul!, roots, isroot_of_unity, promote_rule,
                AbstractAlgebra

export abelian_closure,
       QabField,
       QabElem,
       set_variable!,
       get_variable

################################################################################
#
#  Types
#
################################################################################

mutable struct QabField <: Nemo.Field # union of cyclotomic fields
  s::String
  fields::Dict{Int, AnticNumberField} # Cache for the cyclotomic fields
end

const _Qab = QabField("ζ", Dict{Int, AnticNumberField}())

mutable struct QabElem <: Nemo.FieldElem
  data::nf_elem                       # Element in cyclotomic field
  c::Int                              # Conductor of field
end

# This is a functor like object G with G(n) = primitive n-th root of unity

mutable struct QabFieldGen
  K::QabField
end

const _QabGen = QabFieldGen(_Qab)

################################################################################
#
#  Creation of the field
#
################################################################################

"""
    abelian_closure(QQ::RationalField)

Return a pair `(K, z)` consisting of the abelian closure `K` of the rationals
and a generator `z` that can be used to construct primitive roots of unity in
`K`.
"""
abelian_closure(::FlintRationalField) = _Qab, _QabGen

"""
    gen(K::QabField)

Return the generator of the abelian closure `K` that can be used to construct
primitive roots of unity.
"""
gen(K::QabField) = _QabGen

"""
    gen(K::QabField, s::String)

Return the generator of the abelian closure `K` that can be used to construct
primitive roots of unity. The string `s` will be used during printing.
"""
function gen(K::QabField, s::String)
  K.s = s
  return gen(K)
end

################################################################################
#
#  Parent and element functions
#
################################################################################

elem_type(::Type{QabField}) = QabElem

elem_type(::QabField) = QabElem

parent_type(::Type{QabElem}) = QabField

parent_type(::QabElem) = QabField

Oscar.parent(::QabElem) = _Qab

################################################################################
#
#  Field access
#
################################################################################

_variable(K::QabField) = K.s

_variable(b::QabElem) = Expr(:call, Symbol(_variable(_Qab)), b.c)

function cyclotomic_field(K::QabField, c::Int)
  if haskey(K.fields, c)
    k = K.fields[c]
    return k, gen(k)
  else
    k, z = CyclotomicField(c, string(K.s, "(", c, ")"), cached = false)
    K.fields[c] = k
    return k, z
  end
end

data(a::QabElem) = a.data

################################################################################
#
#  Creation of elements
#
################################################################################

# This function finds a primitive root of unity in our field, note this is
# not always e^(2*pi*i)/n

function root_of_unity(K::QabField, n::Int)
  if n % 2 == 0 && n % 4 != 0
    c = div(n, 2)
  else
    c = n
  end
  K, z = cyclotomic_field(K, c)
  if c == n
    return QabElem(z, c)
  else
    return QabElem(-z, c)
  end
end

function root_of_unity2(K::QabField, c::Int)
  # This function returns the primitive root of unity e^(2*pi*i/n)
  K, z = cyclotomic_field(K, c)
  return QabElem(z, c)
end

(z::QabFieldGen)(n::Int) = root_of_unity(z.K, n)

one(K::QabField) = _Qab(1)

one(a::QabElem) = one(parent(a))

function isone(a::QabElem)
  return isone(data(a))
end

function iszero(a::QabElem)
  return iszero(data(a))
end

zero(K::QabField) = _Qab(0)

zero(a::QabElem) = zero(parent(a))

function (K::QabField)(a::Union{fmpz, fmpq, Integer, Rational})
  return a*root_of_unity(K, 1)
end

function (K::QabField)(a::QabElem)
  return a
end

(K::QabField)() = zero(K)

################################################################################
#
#  String I/O
#
################################################################################

function Base.show(io::IO, a::QabField)
  print(io, "Abelian closure of Q")
end

function Base.show(io::IO, a::QabFieldGen)
  print(io, "Generator of abelian closure of Q")
end

"""
    set_variable!(K::QabField, s::String)

Change the printing of the primitive n-th root of the abelian closure of the
rationals to `s(n)`, where `s` is the supplied string.
"""
function set_variable!(K::QabField, s::String)
  ss = K.s
  K.s = s
  return ss
end

"""
    get_variable(K::QabField)

Return the string used to print the primitive n-th root of the abelian closure
of the rationals.
"""
get_variable(K::QabField) = _variable(K)

function AbstractAlgebra.expressify(b::QabElem; context = nothing)
  a = data(b)
  return AbstractAlgebra.expressify(parent(parent(a).pol)(a), _variable(b), context = context)
end

Oscar.@enable_all_show_via_expressify QabElem

################################################################################
#
#  Singular ring
#
################################################################################

function Oscar.singular_ring(F::QabField)
  return Singular.CoefficientRing(F)
end

################################################################################
#
#  Coercion between cyclotomic fields
#
################################################################################

function isconductor(n::Int)
  if isodd(n)
    return true
  end
  return n % 4 == 0
end

function coerce_up(K::AnticNumberField, n::Int, a::QabElem)
  d = div(n, a.c)
  @assert n % a.c == 0
  #z_n^(d) = z_a
  R = parent(parent(data(a)).pol)
  return QabElem(evaluate(R(data(a)), gen(K)^d), n)
end

function coerce_down(K::AnticNumberField, n::Int, a::QabElem)
  throw(Hecke.NotImplemented())
end

function make_compatible(a::QabElem, b::QabElem)
  if a.c == b.c
    return a,b
  end
  d = lcm(a.c, b.c)
  K, = cyclotomic_field(parent(a), d)
  return coerce_up(K, d, a), coerce_up(K, d, b)
end

################################################################################
#
#  Ring interface functions
#
################################################################################

isunit(a::QabElem) = !iszero(a)

################################################################################
#
#  Arithmetic
#
################################################################################

function +(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(data(a) + data(b), a.c)
end

function *(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(data(a) * data(b), a.c)
end

function -(a::QabElem)
  return QabElem(-data(a), a.c)
end

function ^(a::QabElem, n::Integer)
  return QabElem(data(a)^n, a.c)
end

function ^(a::QabElem, n::fmpz)
  return a^Int(n)
end

function -(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data-b.data, a.c)
end

function //(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data//b.data, a.c)
end

function div(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(a.data//b.data, a.c)
end

function divexact(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return QabElem(divexact(a.data, b.data), a.c)
end

function inv(a::QabElem)
  return QabElem(inv(data(a)), a.c)
end

################################################################################
#
#  Unsafe operations
#
################################################################################

function addeq!(c::QabElem, a::QabElem)
  _c, _a = make_compatible(c, a)
  addeq!(_c.data, _a.data)
  return _c
  #c.data = _c.data
  #c.c = _c.c
  #return c
end

function neg!(a::QabElem)
  mul!(a.data,a.data,-1)
  return a
end

function mul!(c::QabElem, a::QabElem, b::QabElem)
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

*(a::fmpz, b::QabElem) = QabElem(b.data*a, b.c)

*(a::fmpq, b::QabElem) = QabElem(b.data*a, b.c)

*(a::Integer, b::QabElem) = QabElem(data(b) * a, b.c)

*(a::Rational, b::QabElem) = QabElem(data(b) * a, b.c)

*(a::QabElem, b::fmpz) = b*a

*(a::QabElem, b::fmpq) = b*a

*(a::QabElem, b::Integer) = b*a

*(a::QabElem, b::Rational) = b*a

+(a::fmpz, b::QabElem) = QabElem(b.data + a, b.c)

+(a::fmpq, b::QabElem) = QabElem(b.data + a, b.c)

+(a::Integer, b::QabElem) = QabElem(data(b) + a, b.c)

+(a::Rational, b::QabElem) = QabElem(data(b) + a, b.c)

+(a::QabElem, b::fmpz) = b + a

+(a::QabElem, b::fmpq) = b + a

+(a::QabElem, b::Integer) = b + a

+(a::QabElem, b::Rational) = b + a

-(a::fmpz, b::QabElem) = QabElem(-(a, data(b)), b.c)

-(a::fmpq, b::QabElem) = QabElem(-(a, data(b)), b.c)

-(a::Integer, b::QabElem) = QabElem(-(a, data(b)), b.c)

-(a::Rational, b::QabElem) = QabElem(-(a, data(b)), b.c)

-(a::QabElem, b::fmpz) = QabElem(-(data(a), b), a.c)

-(a::QabElem, b::fmpq) = QabElem(-(data(a), b), a.c)

-(a::QabElem, b::Integer) = QabElem(-(data(a), b), a.c)

-(a::QabElem, b::Rational) = QabElem(-(data(a), b), a.c)

//(a::QabElem, b::fmpz) = QabElem(data(a)//b, a.c)

//(a::QabElem, b::fmpq) = QabElem(data(a)//b, a.c)

//(a::QabElem, b::Integer) = QabElem(data(a)//b, a.c)

//(a::QabElem, b::Rational) = QabElem(data(a)//b, a.c)

################################################################################
#
#  Comparison
#
################################################################################

function ==(a::QabElem, b::QabElem)
  a, b = make_compatible(a, b)
  return a.data == b.data
end

function ==(a::QabElem, b::Union{fmpz, fmpq, Integer, Rational})
  return data(a) == b
end

function ==(a::Union{fmpz, fmpq, Integer, Rational}, b::QabElem)
  return b == a
end

################################################################################
#
#  Copy and deepcopy
#
################################################################################

function Base.copy(a::QabElem)
  return QabElem(data(a), a.c)
end

function Base.deepcopy_internal(a::QabElem, dict::IdDict)
  return QabElem(deepcopy_internal(data(a), dict), a.c)
end

#Oscar.isnegative(::QabElem) = false

################################################################################
#
#  Promotion rules
#
################################################################################

AbstractAlgebra.promote_rule(::Type{QabElem}, ::Type{Int}) = QabElem

AbstractAlgebra.promote_rule(::Type{QabElem}, ::Type{fmpz}) = QabElem

AbstractAlgebra.promote_rule(::Type{QabElem}, ::Type{fmpq}) = QabElem

#Oscar.promote_rule(::Type{QabElem}, ::Type{fmpq_poly}) = QabElem

###############################################################################
#
#   Functions for computing roots
#
###############################################################################

function Oscar.root(a::QabElem, n::Int)
  Hecke.@req isroot_of_unity(a) "Element must be a root of unity"
  o = Oscar.order(a)
  l = o*n
  mu = root_of_unity2(parent(a), Int(l))
  return mu
end

function Oscar.roots(a::QabElem, n::Int)
  #Assuming that a is a root of unity, finds all its n-th roots.
  #all roots in a probably smaller field than with the function allrootNew
  #(using root_of_unity -> construct field Q(z_5) when needing a 10th root of unity,
  #root_of_unity constructs the field Q(z_10))
  o = order(a)
  l = o*n
  mu = root_of_unity(parent(a), Int(l))
  A = QabElem[]
  if l==1 && mu==a
    push!(A, mu)
  end
  for k = 0:(l-1)
    el = mu^k
    if el^n == a
      push!(A, el)
    end
  end
  return A
end

function isroot_of_unity(a::QabElem)
  return istorsion_unit(a.data)
  #=
  b = a^a.c
  return b.data == 1 || b.data == -1
  =#
end

function Oscar.order(a::QabElem)
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
#   Galois automorphisms of Qab
#
################################################################################

# The Galois automorphisms of the $n$-th cyclotomic field are the maps
# defined by $\zeta_n \mapsto \zeta_n^k$, for $1 \leq k < n$,
# with $\gcd(n, k) = 1$.
# Thus we can define automorphisms $\sigma_k$ of Qab as follows.
# For each prime power $q$, $\zeta_q$ is mapped to $\zeta_q^k$ if
# $k$ and $q$ are coprime, and to $\zeta_q$ otherwise.
#
# The action of such a map $\sigma_k$ on the $n$-th cyclotomic field can be
# described by $\sigma_l$, with $l$ coprime to $n$:
# Write $n = n_0 n_1$ where $\gcd(n_0, n_1) = 1$ and $n_1$ is maximal
# with $\gcd(k, n_1) = 1$, and choose $a, b$ with $1 = a n_0 + b n_1$.
# Then $l = k a n_0 + b n_1$ is coprime to $n$ and has the properties
# $l \equiv 1 \pmod{n_0}$ and $l \equiv k \pmod{n_1}$.

mutable struct QabAutomorphism
  exp::Int
end

Oscar.hom(K::QabField, L::QabField, k::Int) = QabAutomorphism(k)

(f::QabAutomorphism)(a::QabElem) = a^f

function ^(val::QabElem, sigma::QabAutomorphism)
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
  return QabElem(F(R(res)), n)
end

###############################################################################
#
#   Elements in quadratic subfields of cyclotomic fields
#
###############################################################################

function generators_galois_group_cyclotomic_field(n::Int)
  res = GAP.Globals.GeneratorsPrimeResidues(GAP.julia_to_gap(n))
  return [QabAutomorphism(k)
          for k in Vector{Int}(GAP.Globals.Flat(res.generators))]
end

"""
    square_root_in_cyclotomic_field(F::QabField, n::Int, N::Int)

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
function square_root_in_cyclotomic_field(F::QabField, n::Int, N::Int)
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

  return cf * QabElem(FF(elm), N)
end

"""
    quadratic_irrationality_info(a::QabModule.QabElem)

Return `(x, y, n)`, where `x`, `y` are of type `fmpq` and `n` is
a squarefree integer, such that `a == x + y sqrt(n)` holds.

(We assume that the underlying primitive `N`-th root of unity that
is used to define `a` is identified with the complex number `exp(2*Pi*i/N)`,
where `i` is the imaginary unit.)
"""
function quadratic_irrationality_info(a::QabElem)
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
    for sigma in galgens
      img = cand^sigma
      if img != a && img != cand
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

using .AbelianClosure

export abelian_closure,
       QabField,
       QabElem,
       set_variable!,
       get_variable

