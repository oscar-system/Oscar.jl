export restrict, generic_fraction, restriction_map, is_identity_map, canonical_isomorphism 
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
    V::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing}
  )
  isempty(V) && return zero(OO(V))
  for i in 1:length(restrictions(f))
    if V === affine_patches(domain(f))[i]
      return restrictions(f)[i]
    end
  end
  issubset(V, domain(f)) || error("the set is not contained in the domain of definition of the function")
  VU = [intersect(V, U) for U in affine_patches(domain(f))]
  g = [OO(VU[i])(f[i]) for i in 1:length(VU)]
  l = write_as_linear_combination(one(OO(V)), OO(V).(lifted_denominator.(g)))
  a = dot(l, OO(V).(lifted_numerator.(g)))
  return a
end

function restrict(
    f::SpecOpenRingElem, 
    V::AbsSpec{<:Ring, <:MPolyLocalizedRing}
  )
  isempty(V) && return zero(OO(V))
  for i in 1:length(restrictions(f))
    if V === affine_patches(domain(f))[i]
      return restrictions(f)[i]
    end
  end
  issubset(V, domain(f)) || error("the set is not contained in the domain of definition of the function")
  VU = [intersect(V, U) for U in affine_patches(domain(f))]
  g = [OO(VU[i])(f[i]) for i in 1:length(VU)]
  J = ideal(OO(V), denominator.(g))
  l = coordinates(one(OO(V)), ideal(OO(V), denominator.(g)))
  a = dot(l, OO(V).(numerator.(g)))
  return a
end

function restrict(
    f::SpecOpenRingElem, 
    V::SpecOpen
  )
  V === domain(parent(f)) && return f
  fres = [restrict(f, V[i]) for i in 1:ngens(V)]
  return SpecOpenRingElem(OO(V), fres, check=false)
end

########################################################################
# Generic fractions                                                    #
########################################################################
@Markdown.doc """
    generic_fraction(a::SpecOpenRingElem, U::SpecOpen)

Given a regular function ``a ∈ 𝒪(U)`` on a Zariski open 
subset ``U ⊂ X`` of an affine scheme ``X``, return a 
fraction ``p/q`` in `Quot(P)` (where ``P`` is the `ambient_ring` of 
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
  X = ambient(U)
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
function ^(a::SpecOpenRingElem, i::fmpz)
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

### TODO: Rethink this. For instance, restrictions can happen both from and to Specs.
function AbstractAlgebra.promote_rule(::Type{T}, ::Type{RET}) where {T<:SpecOpenRingElem, RET<:RingElem} 
  return T
end

########################################################################
# Additional methods for compatibility and coherence                   #
########################################################################
function (R::MPolyQuo)(a::RingElem, b::RingElem; check::Bool=true)
  return R(a)*inv(R(b))
end

function (R::MPolyRing)(a::RingElem, b::RingElem; check::Bool=true)
  return R(a)*inv(R(b))
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
    h::MPolyElem; 
    check::Bool=true
  )
  Y = ambient(U)

  # handle the shortcut 
  if any(x->(x===X), affine_patches(U))
    i = findfirst(x->(x === X), affine_patches(U))
    function mymap(f::SpecOpenRingElem)
      return f[i]
    end
    return MapFromFunc(mymap, OO(U), OO(X))
  end

  # do the checks
  if check
    X == hypersurface_complement(Y, h) || error("$X is not the hypersurface complement of $h in the ambient variety of $U")
    issubset(X, U) || error("$X is not a subset of $U")
  end

  # first find some basic relation hᵏ= ∑ᵢ aᵢ⋅dᵢ
  d = gens(U)
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
  S, t = PolynomialRing(W, ["t$i" for i in 1:r])
  ta = sum([t*a for (t, a) in zip(t, a)])
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
    for (b, m) in zip(coefficients(dirty), monomials(dirty))
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
      for (b, m) in zip(coefficients(dirty), monomials(dirty))
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
    g = [W(p, q, check=false) for (p, q, dk, k) in sep]
    dk = [dk for (p, q, dk, k) in sep]
    return OO(X)(sum([a*b for (a, b) in zip(g, c)]), check=false)*OO(X)(1//poh^m, check=false)
  end
  return Hecke.MapFromFunc(mysecondmap, OO(U), OO(X))
end

# Automatically find a hypersurface equation h such that X = D(h) in 
# the ambient scheme Y of U. 
function restriction_map(U::SpecOpen{<:AbsSpec{<:Ring, <:AbsLocalizedRing}},
    X::AbsSpec{<:Ring, <:AbsLocalizedRing}; 
    check::Bool=true
  )
  Y = ambient(U)
  R = ambient_ring(Y)
  R == ambient_ring(X) || error("`ambient_ring`s of the schemes not compatible")
  if check
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
function restriction_map(U::SpecOpen{<:AbsSpec{<:Ring, <:MPolyQuo}},
    X::AbsSpec{<:Ring, <:AbsLocalizedRing};
    check::Bool=true
  )
  Y = ambient(U)
  R = ambient_ring(Y)
  R == ambient_ring(X) || error("rings not compatible")
  if check
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
  Y = ambient(U)
  R = ambient_ring(Y)
  R == ambient_ring(X) || error("rings not compatible")
  if check
    issubset(X, Y) || error("$X is not contained in the ambient scheme of $U")
    issubset(X, U) || error("$X is not a subset of $U")
  end
  h = prod(denominators(inverted_set(OO(X))))
  return restriction_map(U, X, h, check=false)
end


# For f = p//q and d this computes a decomposition p', q', d^k, k 
# such that f = p'//(q'⋅d^k) and q' and d have no common factors. 
function pull_from_denominator(f::MPolyQuoLocalizedRingElem, d::MPolyElem)
  p = lifted_numerator(f)
  q = lifted_denominator(f)
  (i, o) = ppio(q, d)
  (k, pod) = Oscar._minimal_power_such_that(d, x->(divides(x, i)[1]))
  b = divexact(pod, i)
  return b*p, o, pod, k
end

function pull_from_denominator(f::MPolyLocalizedRingElem, d::MPolyElem)
  p = numerator(f)
  q = denominator(f)
  (i, o) = ppio(q, d)
  (k, pod) = Oscar._minimal_power_such_that(d, x->(divides(x, i)[1]))
  b = divexact(pod, i)
  return b*p, o, pod, k
end

function restriction_map(X::Spec, U::SpecOpen; check::Bool=true)
  Y = ambient(U)
  if check
    all(V->issubset(V, X), affine_patches(U)) || error("$U is not a subset of $X")
  end
  function mymap(f::MPolyQuoLocalizedRingElem)
    return SpecOpenRingElem(OO(U), [OO(V)(f) for V in affine_patches(U)])
  end
  return Hecke.MapFromFunc(mymap, OO(X), OO(U))
end

function restriction_map(U::SpecOpen, V::SpecOpen; check::Bool=true)
  if check
    issubset(V, U) || error("$V is not a subset of $U")
  end

  if U === V
    function mymap(f::SpecOpenRingElem)
      return f
    end
    return Hecke.MapFromFunc(mymap, OO(U), OO(V))
  end

  if ambient(U) === ambient(V)
    g = [restriction_map(U, W, d, check=false) for (W, d) in zip(affine_patches(V), gens(V))]
    function mysecondmap(f::SpecOpenRingElem)
      return SpecOpenRingElem(OO(V), [h(f) for h in g], check=false)
    end
    return Hecke.MapFromFunc(mysecondmap, OO(U), OO(V))
  end
  
  g = [restriction_map(U, W, check=false) for W in affine_patches(V)]
  function mythirdmap(f::SpecOpenRingElem)
    return SpecOpenRingElem(OO(V), [g(f) for g in g], check=false)
  end
  return Hecke.MapFromFunc(mythirdmap, OO(U), OO(V))
end

########################################################################
# Maps of SpecOpenRings                                                #
########################################################################
function is_identity_map(f::Hecke.Map{DomType, CodType}) where {DomType<:SpecOpenRing, CodType<:SpecOpenRing}
  domain(f) === codomain(f) || return false
  R = ambient_ring(scheme(domain(f)))
  return all(x->(domain(f)(x) == f(domain(f)(x))), gens(R))
end

function canonical_isomorphism(S::SpecOpenRing, T::SpecOpenRing; check::Bool=true)
  X = scheme(S)
  Y = scheme(T)
  R = ambient_ring(X)
  R == ambient_ring(Y) || error("rings can not be canonically compared")
  if check
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
  return Hecke.MapFromFunc(mymap, myinvmap, S, T)
end

