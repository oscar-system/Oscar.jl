########################################################################
# Methods for SpecOpenRing                                             #
########################################################################

function ==(R::SpecOpenRing, S::SpecOpenRing)
  return scheme(R)==scheme(S) && domain(R)==domain(S)
end

########################################################################
# Methods for SpecOpenRingElem                                         #
########################################################################

########################################################################
# Restrictions of regular functions                                    #
########################################################################
function restrict(
    f::SpecOpenRingElem, 
    V::AbsSpec{<:Ring, <:MPolyQuoLocRing};
    check::Bool=true
  )
  check && isempty(V) && return zero(OO(V))
  for i in 1:length(restrictions(f))
    if V === affine_patches(domain(f))[i]
      return restrictions(f)[i]
    end
  end
  check && (issubset(V, domain(f)) || error("the set is not contained in the domain of definition of the function"))
  VU = [intersect(V, U) for U in affine_patches(domain(f))]
  g = [OO(VU[i])(f[i]) for i in 1:length(VU)]
  #l = write_as_linear_combination(one(OO(V)), OO(V).(lifted_denominator.(g)))
  l = coordinates(one(OO(V)), ideal(OO(V), lifted_denominator.(g)))
  a = dot(l, OO(V).(lifted_numerator.(g)))
  return a
end

function restrict(
    f::SpecOpenRingElem, 
    V::AbsSpec{<:Ring, <:MPolyLocRing};
    check::Bool=true
  )
  check && isempty(V) && return zero(OO(V))
  for i in 1:length(restrictions(f))
    if V === affine_patches(domain(f))[i]
      return restrictions(f)[i]
    end
  end
  check && (issubset(V, domain(f)) || error("the set is not contained in the domain of definition of the function"))
  VU = [intersect(V, U) for U in affine_patches(domain(f))]
  g = [OO(VU[i])(f[i]) for i in 1:length(VU)]
  J = ideal(OO(V), denominator.(g))
  l = coordinates(one(OO(V)), ideal(OO(V), denominator.(g)))
  a = dot(l, OO(V).(numerator.(g)))
  return a
end

function restrict(
    f::SpecOpenRingElem, 
    V::SpecOpen; 
    check::Bool=true # Only for compatibility of the call signature
  )
  V === domain(parent(f)) && return f
  fres = [restrict(f, V[i]) for i in 1:ngens(V)]
  return SpecOpenRingElem(OO(V), fres, check=false)
end

########################################################################
# Generic fractions                                                    #
########################################################################
@doc raw"""
    generic_fraction(a::SpecOpenRingElem, U::SpecOpen)

Given a regular function ``a ∈ 𝒪(U)`` on a Zariski open 
subset ``U ⊂ X`` of an affine scheme ``X``, return a 
fraction ``p/q`` in `Quot(P)` (where ``P`` is the `ambient_coordinate_ring` of
the `ambient` scheme ``X`` of ``U``) which represents ``a``
in the sense that the maximal extension of its restriction
to ``U`` returns ``a``.

**Note:** The seemingly superfluous argument ``U`` is needed 
to have a coherent syntax with the method for regular functions 
``a`` on `PrincipalOpenSubset`s. There, the element ``a`` does 
not know about its scheme, so it has to be passed as an extra argument.
"""
function generic_fraction(a::SpecOpenRingElem, U::SpecOpen)
  U === domain(a) || error("domains are not compatible")
  X = ambient_scheme(U)
  d = find_non_zero_divisor(U)
  W = hypersurface_complement(X, d)
  b = restrict(a, W)
  return lifted_numerator(b)//lifted_denominator(b)
end

########################################################################
# Arithmetic                                                           #
########################################################################
function +(a::T, b::T) where {T<:SpecOpenRingElem}
  parent(a) === parent(b) || return a + (parent(a)(b))
  return SpecOpenRingElem(parent(a), [a[i] + b[i] for i in 1:length(restrictions(a))], check=false)
end

function -(a::T, b::T) where {T<:SpecOpenRingElem}
  parent(a) === parent(b) || return a - (parent(a)(b))
  return SpecOpenRingElem(parent(a), [a[i] - b[i] for i in 1:length(restrictions(a))], check=false)
end

function -(a::T) where {T<:SpecOpenRingElem}
  return SpecOpenRingElem(parent(a), [-a[i] for i in 1:length(restrictions(a))], check=false)
end

function *(a::T, b::T) where {T<:SpecOpenRingElem}
  parent(a) === parent(b) || return a * (parent(a)(b))
  return SpecOpenRingElem(parent(a), [a[i] * b[i] for i in 1:length(restrictions(a))], check=false)
end

#function *(a::RingElem, b::T) where {T<:SpecOpenRingElem}
#  return b*(parent(b)(a))
#end

function *(a::Integer, b::T) where {T<:SpecOpenRingElem}
  return b*(parent(b)(a))
end

#function *(b::T, a::RingElem) where {T<:SpecOpenRingElem}
#  return a*b
#end

function ==(a::T, b::T) where {T<:SpecOpenRingElem}
  parent(a) === parent(b) || return a == (parent(a)(b))
  for i in 1:length(restrictions(a))
    a[i] == b[i] || return false
  end
  return true
end

function ^(a::SpecOpenRingElem, i::Int64)
  return SpecOpenRingElem(parent(a), [a[k]^i for k in 1:length(restrictions(a))])
end
function ^(a::SpecOpenRingElem, i::Integer)
  return SpecOpenRingElem(parent(a), [a[k]^i for k in 1:length(restrictions(a))])
end
function ^(a::SpecOpenRingElem, i::ZZRingElem)
  return SpecOpenRingElem(parent(a), [a[k]^i for k in 1:length(restrictions(a))])
end

function divexact(a::T, b::T; check::Bool=false) where {T<:SpecOpenRingElem} 
  parent(a) === parent(b) || return divexact(a, (parent(a)(b)))
  return SpecOpenRingElem(parent(a), [divexact(a[i], b[i]) for i in 1:length(restrictions(a))])
end

function is_unit(a::SpecOpenRingElem) 
  return all(x->is_unit(x), restrictions(a))
end

inv(a::SpecOpenRingElem) = SpecOpenRingElem(parent(a), [inv(f) for f in restrictions(a)], check=false)

########################################################################
# Promote rules for the ring interface                                 #
########################################################################
AbstractAlgebra.promote_rule(::Type{T}, ::Type{RET}) where {T<:SpecOpenRingElem, RET<:Integer} = T
AbstractAlgebra.promote_rule(::Type{RET}, ::Type{T}) where {T<:SpecOpenRingElem, RET<:Integer} = T

### Promote rules from and to other rings can not be made in a coherent 
# way depending only on the types. One problem is that both, restricting 
# from and to an affine scheme are valid operations, depending on the 
# specific geometric configuration. Hence, we rely on the user to perform 
# every coercion manually. 

# Additional promotions to make (graded) polynomial rings and their quotients work over SpecOpenRings.
*(a::S, b::T) where {S<:SpecOpenRingElem, T<:MPolyRingElem{S}} = parent(b)(a)*b
*(a::S, b::T) where {S<:SpecOpenRingElem, T<:MPolyDecRingElem{S}} = parent(b)(a)*b
*(a::S, b::T) where {S<:SpecOpenRingElem, T<:MPolyQuoRingElem{MPolyDecRingElem{S}}} = parent(b)(a)*b
*(a::S, b::T) where {S<:SpecOpenRingElem, T<:MPolyQuoRingElem{MPolyRingElem{S}}} = parent(b)(a)*b

*(b::T, a::S) where {S<:SpecOpenRingElem, T<:MPolyDecRingElem{S}} = parent(b)(a)*b
*(b::T, a::S) where {S<:SpecOpenRingElem, T<:MPolyQuoRingElem{MPolyDecRingElem{S}}} = parent(b)(a)*b
*(b::T, a::S) where {S<:SpecOpenRingElem, T<:MPolyQuoRingElem{MPolyRingElem{S}}} = parent(b)(a)*b

########################################################################
# Additional methods for compatibility and coherence                   #
########################################################################
function _cast_fraction(R::MPolyQuoRing, a::RingElem, b::RingElem; check::Bool=true)
  return R(a)*inv(R(b))
end

function _cast_fraction(R::MPolyRing, a::RingElem, b::RingElem; check::Bool=true)
  return R(a)*inv(R(b))
end

function _cast_fraction(R::Union{<:MPolyLocRing, <:MPolyQuoLocRing}, a, b; check::Bool=true)
  return R(a, b; check)
end

########################################################################
# Restrictions of regular functions from SpecOpens to AbsSpecs         #
# and vice versa                                                       #
########################################################################

# For an open subset U ⊂ Y of an affine scheme Y and a hypersurface
# complement X = D(h) ⊂ Y with X ⊂ U this returns the restriction
# map ρ : 𝒪(U) → 𝒪(X)
function restriction_map(
    U::SpecOpen, 
    X::AbsSpec{<:Ring, <:AbsLocalizedRing}, 
    h::MPolyRingElem; 
    check::Bool=true
  )
  Y = ambient_scheme(U)

  # handle the shortcut 
  if any(x->(x===X), affine_patches(U))
    i = findfirst(x->(x === X), affine_patches(U))
    function mymap(f::SpecOpenRingElem)
      return f[i]
    end
    return MapFromFunc(OO(U), OO(X), mymap)
  end

  # do the checks
  @check begin
    X == hypersurface_complement(Y, h) || error("$X is not the hypersurface complement of $h in the ambient variety of $U")
    issubset(X, U) || error("$X is not a subset of $U")
  end

  # first find some basic relation hᵏ= ∑ᵢ aᵢ⋅dᵢ
  d = complement_equations(U)
  I = complement_ideal(U)
  # _minimal_power_such_that(P, h) returns a tuple (k, h^k) with 
  # k the minimal exponent such that the property P(h^k) returns `true`.
  (k, poh) = Oscar._minimal_power_such_that(h, x->(base_ring(I)(x) in I))
  a = coordinates(base_ring(I)(poh), I)
  r = length(d)

  # the local representatives of the input f will be of the form gᵢ⋅1//dᵢˢ⁽ⁱ⁾
  # with gᵢ∈ 𝒪(Y). For higher powers s(i) > 1 we need other coefficients 
  # cᵢ for the relation 
  #
  #   hˡ = ∑ᵢ cᵢ⋅dˢ⁽ⁱ⁾
  #
  # for some power hˡ. To this end, we set up a polynomial ring 𝒪(Y)[t₁,…,tᵣ]
  # and take powers of the element ∑ᵢaᵢ⋅tᵢ with the coefficients aᵢ of the basic 
  # relation. Eventually, all terms appearing in that expression will have 
  # monomials t₁ᵉ⁽¹⁾⋅…⋅tᵣᵉ⁽ʳ⁾ with some e(i) ≥ s(i). Substituting and grouping 
  # the terms accordingly, we derive the desired expressions for the cᵢ's.
  #W = localized_ring(OO(Y))
  W = OO(Y)
  S, t = polynomial_ring(W, ["t$i" for i in 1:r])
  ta = length(a) == 0 ? zero(S) : sum([t*a for (t, a) in zip(t, a)])
  function mysecondmap(f::SpecOpenRingElem)
    sep = [pull_from_denominator(f[i], d[i]) for i in 1:r]
    # the following takes care of oddities from zero divisors.
    for i in 1:r-1
      for j in i+1:r
        while !iszero(OO(Y)(sep[i][1]*sep[j][2]*sep[j][3] - sep[j][1]*sep[i][2]*sep[i][3]))
          sep[i] = (sep[i][1]*d[i], sep[i][2], sep[i][3]*d[i], sep[i][4] + 1)
          sep[j] = (sep[j][1]*d[j], sep[j][2], sep[j][3]*d[j], sep[j][4] + 1)
        end
      end
    end

    k = [k for (p, q, dk, k) in sep]
    c = [zero(W) for i in 1:r]
    dirty = one(S)
    m = 0
    # one extra round to catch the degenerate case where no powers are needed
    cleaned = zero(dirty)
    for (b, m) in zip(AbstractAlgebra.coefficients(dirty), AbstractAlgebra.monomials(dirty))
      for i in 1:r
        if exponent(m, 1, i) == k[i]
          c[i] = c[i] + b*evaluate(m, [(j == i ? one(W) : W(d[j])) for j in 1:r])
          cleaned = cleaned + b*m
          break
        end
      end
    end
    dirty = dirty - cleaned

    while !iszero(dirty)
      m = m + 1
      c = (x->poh*x).(c)
      dirty = dirty*ta
      cleaned = zero(dirty)
      for (b, m) in zip(AbstractAlgebra.coefficients(dirty), AbstractAlgebra.monomials(dirty))
        for i in 1:r
          if exponent(m, 1, i) == k[i]
            c[i] = c[i] + b*evaluate(m, [(j == i ? one(W) : W(d[j])) for j in 1:r])
            cleaned = cleaned + b*m
            break
          end
        end
      end
      dirty = dirty - cleaned
    end
    g = [_cast_fraction(W, p, q, check=false) for (p, q, dk, k) in sep]
    dk = [dk for (p, q, dk, k) in sep]
    return OO(X)(sum([a*b for (a, b) in zip(g, c)]), check=false)*OO(X)(1//poh^m, check=false)
  end
  return MapFromFunc(OO(U), OO(X), mysecondmap)
end

# Automatically find a hypersurface equation h such that X = D(h) in 
# the ambient scheme Y of U. 
function restriction_map(U::SpecOpen{<:AbsSpec{<:Ring, <:AbsLocalizedRing}},
    X::AbsSpec{<:Ring, <:AbsLocalizedRing}; 
    check::Bool=true
  )
  Y = ambient_scheme(U)
  R = ambient_coordinate_ring(Y)
  R == ambient_coordinate_ring(X) || error("`ambient_coordinate_ring`s of the schemes not compatible")
  @check begin
    issubset(X, Y) || error("$X is not contained in the ambient scheme of $U")
    issubset(X, U) || error("$X is not a subset of $U")
  end
  L = localized_ring(OO(X))
  D = denominators(inverted_set(L))
  p = prod(denominators(inverted_set(OO(Y))))
  h = one(R)
  for d in D
    (i, o) = ppio(d, p)
    h = h*o
  end
  return restriction_map(U, X, h, check=false)
end

# Automatically find a hypersurface equation h such that X = D(h) in
# the ambient scheme Y of U.
function restriction_map(U::SpecOpen{<:AbsSpec{<:Ring, <:MPolyQuoRing}},
    X::AbsSpec{<:Ring, <:AbsLocalizedRing};
    check::Bool=true
  )
  Y = ambient_scheme(U)
  R = ambient_coordinate_ring(Y)
  R == ambient_coordinate_ring(X) || error("rings not compatible")
  @check begin
    issubset(X, Y) || error("$X is not contained in the ambient scheme of $U")
    issubset(X, U) || error("$X is not a subset of $U")
  end
  h = prod(denominators(inverted_set(OO(X))))
  return restriction_map(U, X, h, check=false)
end

function restriction_map(U::SpecOpen{<:AbsSpec{<:Ring, <:MPolyRing}},
    X::AbsSpec{<:Ring, <:AbsLocalizedRing};
    check::Bool=true
  )
  Y = ambient_scheme(U)
  R = ambient_coordinate_ring(Y)
  R == ambient_coordinate_ring(X) || error("rings not compatible")
  @check begin
    issubset(X, Y) || error("$X is not contained in the ambient scheme of $U")
    issubset(X, U) || error("$X is not a subset of $U")
  end
  h = prod(denominators(inverted_set(OO(X))))
  return restriction_map(U, X, h, check=false)
end


# For f = p//q and d this computes a decomposition p', q', d^k, k 
# such that f = p'//(q'⋅d^k) and q' and d have no common factors. 
function pull_from_denominator(f::MPolyQuoLocRingElem, d::MPolyRingElem)
  p = lifted_numerator(f)
  q = lifted_denominator(f)
  (i, o) = ppio(q, d)
  (k, pod) = Oscar._minimal_power_such_that(d, x->(divides(x, i)[1]))
  b = divexact(pod, i)
  return b*p, o, pod, k
end

function pull_from_denominator(f::MPolyLocRingElem, d::MPolyRingElem)
  p = numerator(f)
  q = denominator(f)
  (i, o) = ppio(q, d)
  (k, pod) = Oscar._minimal_power_such_that(d, x->(divides(x, i)[1]))
  b = divexact(pod, i)
  return b*p, o, pod, k
end

function restriction_map(X::Spec, U::SpecOpen; check::Bool=true)
  Y = ambient_scheme(U)
  @check all(V->issubset(V, X), affine_patches(U)) "$U is not a subset of $X"
  function mymap(f::MPolyQuoLocRingElem)
    return SpecOpenRingElem(OO(U), [OO(V)(f) for V in affine_patches(U)])
  end
  return MapFromFunc(OO(X), OO(U), mymap)
end

function restriction_map(U::SpecOpen, V::SpecOpen; check::Bool=true)
  @check issubset(V, U) "$V is not a subset of $U"

  if U === V
    function mymap(f::SpecOpenRingElem)
      return f
    end
    return MapFromFunc(OO(U), OO(V), mymap)
  end

  if ambient_scheme(U) === ambient_scheme(V)
    g = [restriction_map(U, W, d, check=false) for (W, d) in zip(affine_patches(V), complement_equations(V))]
    function mysecondmap(f::SpecOpenRingElem)
      return SpecOpenRingElem(OO(V), [h(f) for h in g], check=false)
    end
    return MapFromFunc(OO(U), OO(V), mysecondmap)
  end
  
  g = [restriction_map(U, W, check=false) for W in affine_patches(V)]
  function mythirdmap(f::SpecOpenRingElem)
    return SpecOpenRingElem(OO(V), [g(f) for g in g], check=false)
  end
  return MapFromFunc(OO(U), OO(V), mythirdmap)
end

########################################################################
# Maps of SpecOpenRings                                                #
########################################################################
function is_identity_map(f::Map{DomType, CodType}) where {DomType<:SpecOpenRing, CodType<:SpecOpenRing}
  domain(f) === codomain(f) || return false
  R = ambient_coordinate_ring(scheme(domain(f)))
  return all(x->(domain(f)(x) == f(domain(f)(x))), gens(R))
end

function canonical_isomorphism(S::SpecOpenRing, T::SpecOpenRing; check::Bool=true)
  X = scheme(S)
  Y = scheme(T)
  R = ambient_coordinate_ring(X)
  R == ambient_coordinate_ring(Y) || error("rings can not be canonically compared")
  @check begin
    (domain(S) == domain(T)) || error("open domains are not isomorphic")
  end

  pb_to_Vs = [restriction_map(domain(S), V) for V in affine_patches(domain(T))]
  pb_to_Us = [restriction_map(domain(T), U) for U in affine_patches(domain(S))]
  function mymap(a::SpecOpenRingElem)
    return SpecOpenRingElem(T, [g(a) for g in pb_to_Vs], check=false)
  end
  function myinvmap(b::SpecOpenRingElem)
    return SpecOpenRingElem(S, [g(b) for g in pb_to_Us], check=false)
  end
  return MapFromFunc(S, T, mymap, myinvmap)
end

# Special override for a case where even ideal membership and ring flattenings 
# are not implemented. 
function simplify(f::MPolyQuoRingElem{<:MPolyDecRingElem{<:SpecOpenRingElem}})
  return f
end

###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, R::SpecOpenRing)
  if get(io, :supercompact, false)
    print(io, "Ring of regular functions")
  else
    io = pretty(io)
    print(io, "Ring of regular functions on ", Lowercase(), domain(R))
  end
end

# Here we just need some details on the domain `U`.
function Base.show(io::IO, ::MIME"text/plain", R::SpecOpenRing)
  io = pretty(io)
  println(io, "Ring of regular functions")
  print(io, Indent(), "on ", Lowercase())
  show(io, domain(R))
  print(io, Dedent())
end

function Base.show(io::IO, a::SpecOpenRingElem)
  if get(io, :supercompact, false)
    print(io, "Reguler function")
  else
    io = pretty(io)
    print(io, "Regular function on ", Lowercase(), domain(parent(a)))
  end
end

# Since regular functions are described on each affine patches of the domain `U`
# on which they are defined, we need to extract details about the affine patches
# on `U` and label them so that one can see how the regular function is defined
# on each patch
function Base.show(io::IO, ::MIME"text/plain", a::SpecOpenRingElem)
  io = pretty(io)
  R = parent(a)
  U = domain(R)
  println(io, "Regular function")
  print(io, Indent(), "on ", Lowercase())
  show(IOContext(io, :show_semi_compact => true), domain(R))
  print(io, Dedent())
  r = restrictions(a)
  ap = affine_patches(a)
  if length(r) > 0
    l = ndigits(length(r))
    println(io)
    print(io, "with restriction")
    length(r) > 1 && print(io, "s")
    print(io, Indent())
    for i in 1:length(r)
      li = ndigits(i)
      println(io)
      print(io, "patch", " "^(l-li+1)*"$(i): ", r[i])
    end
    print(io, Dedent())
  end
end
  
