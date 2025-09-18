################################################################################
# Printing and IO
################################################################################

function Base.show(io::IO, ::MIME"text/plain", P::AbsProjectiveScheme{<:Any, <:MPolyQuoRing})
  io = pretty(io)
  println(io, "Projective scheme")
  println(io, Indent(), "over ", Lowercase(), base_ring(P))
  print(io, Dedent(), "defined by ", Lowercase(), defining_ideal(P))
end

function Base.show(io::IO, P::AbsProjectiveScheme{<:Any, <:MPolyQuoRing})
  if is_terse(io)
    print(io, "Projective scheme")
  elseif get_attribute(P, :is_empty, false)
    io = pretty(io)
    print(io, "Empty projective scheme over ")
    K = base_ring(P)
    print(terse(io), Lowercase(), K)
  else
    io = pretty(io)
    print(io, "Projective scheme in ")
    print(terse(io), Lowercase(), ambient_space(P), " over ", Lowercase(), base_ring(P))
  end
end

# Projective space
function Base.show(io::IO, ::MIME"text/plain", P::AbsProjectiveScheme{<:Any, <:MPolyDecRing})
  io = pretty(io)
  println(io, "Projective space of dimension $(relative_ambient_dimension(P))")
  print(io, Indent(), "over ")
  println(io, Lowercase(), base_ring(P))
  print(io, Dedent(), "with homogeneous coordinate")
  length(homogeneous_coordinates(P)) != 1 && print(io, "s")
  print(io, " [")
  print(io, join(homogeneous_coordinates(P), ", "), "]")
end

function Base.show(io::IO, P::AbsProjectiveScheme{<:Any, <:MPolyDecRing})
  io = pretty(io)
  if is_terse(io)
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, LowercaseOff(), "ℙ$(ltx["\\^$(relative_ambient_dimension(P))"])")
    else
      print(io, LowercaseOff(), "IP^$(relative_ambient_dimension(P))")
    end
  elseif get_attribute(P, :is_empty, false)
    print(io, "Empty projective space over ")
    K = base_ring(P)
    print(terse(io), Lowercase(), K)
  else
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, LowercaseOff(), "ℙ")
      n = relative_ambient_dimension(P)
      for d in reverse(digits(n))
        print(io, ltx["\\^$d"])
      end
      print(io, " over ")
      print(terse(io), Lowercase(), base_ring(P))
    else
      print(io, "Projective $(relative_ambient_dimension(P))-space over ")
      print(terse(io), Lowercase(), base_ring(P))
      c = homogeneous_coordinates(P)
      print(io, " with coordinate")
      length(c) != 1 && print(io, "s")
      print(io, " [")
      print(io, join(c, ", "), "]")
    end
  end
end


@doc raw"""
    dehomogenization_map(X::AbsProjectiveScheme, U::AbsAffineScheme)

Return the restriction morphism from the graded coordinate ring of ``X`` to `𝒪(U)`.

# Examples
```jldoctest
julia> P = projective_space(QQ, ["x0", "x1", "x2"])
Projective space of dimension 2
  over rational field
with homogeneous coordinates [x0, x1, x2]

julia> X = covered_scheme(P);

julia> U = first(affine_charts(X))
Spectrum
  of multivariate polynomial ring in 2 variables (x1//x0), (x2//x0)
    over rational field

julia> phi = dehomogenization_map(P, U);

julia> S = homogeneous_coordinate_ring(P);

julia> phi(S[2])
(x1//x0)

```
"""
function dehomogenization_map(X::AbsProjectiveScheme, U::AbsAffineScheme)
  cache = _dehomogenization_cache(X)
  if haskey(cache, U)
    return cache[U]
  end
  Y = covered_scheme(X)
  C = default_covering(Y)
  inc, _ = _find_chart(U, C)
  V = codomain(inc)
  @assert V in keys(_dehomogenization_cache(X)) "dehomogenization map not found"
  return compose(dehomogenization_map(X, V), OO(Y)(V, U))
end

function _dehomogenization_map(X::AbsProjectiveScheme, U::AbsAffineScheme, i::Int)
  S = homogeneous_coordinate_ring(X)
  s = vcat(gens(OO(U))[1:i-1], [one(OO(U))], gens(OO(U))[i:relative_ambient_dimension(X)])
  phi = hom(S, OO(U), s, check=false)
  return phi
end

function _dehomogenization_map(
    X::AbsProjectiveScheme{CRT},
    U::AbsAffineScheme,
    i::Int
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  S = homogeneous_coordinate_ring(X)
  R = base_ring(X)
  r = relative_ambient_dimension(X)
  p = hom(R, OO(U), gens(OO(U))[r+1:end], check=false)
  s = vcat(gens(OO(U))[1:i-1], [one(OO(U))], gens(OO(U))[i:relative_ambient_dimension(X)])
  phi = hom(S, OO(U), p, s, check=false)
  return phi
end

#=
function dehomogenization_map(
    X::AbsProjectiveScheme{CRT},
    U::AbsAffineScheme
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  return dehomogenization_map(X, X[U][2]-1)
end
=#

#=
@doc raw"""
    dehomogenization_map(X::AbsProjectiveScheme, i::Int)

Return the restriction morphism from the graded coordinate ring of ``X`` to `𝒪(Uᵢ)`.
Where `Uᵢ` is the `i`-th affine chart of `X`.
"""
function dehomogenization_map(X::AbsProjectiveScheme, i::Int)
  i in 0:relative_ambient_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = default_covering(covered_scheme(X))
  U = C[i+1]
  cache = _dehomogenization_cache(X)
  if haskey(cache, U)
    return cache[U]
  end
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  phi = hom(S, OO(U), s, check=false)
  cache[U] = phi
  return phi
end
=#


@doc raw"""
    homogenization_map(P::AbsProjectiveScheme, U::AbsAffineScheme)

Given an affine chart ``U ⊂ P`` of an `AbsProjectiveScheme`
``P``, return a method ``h`` for the homogenization of elements
``a ∈ 𝒪(U)``.

This means that ``h(a)`` returns a pair ``(p, q)`` representing a fraction
``p/q ∈ S`` of the `ambient_coordinate_ring` of ``P`` such
that ``a`` is the dehomogenization of ``p/q``.

**Note:** For the time being, this only works for affine
charts which are of the standard form ``sᵢ ≠ 0`` for ``sᵢ∈ S``
one of the homogeneous coordinates of ``P``.

**Note:** Since this map returns representatives only, it
is not a mathematical morphism and, hence, in particular
not an instance of `Map`.

# Examples
```jldoctest
julia> A, _ = QQ[:u, :v];

julia> P = projective_space(A, ["x0", "x1", "x2"])
Projective space of dimension 2
  over multivariate polynomial ring in 2 variables over QQ
with homogeneous coordinates [x0, x1, x2]

julia> X = covered_scheme(P)
Scheme
  over rational field
with default covering
  described by patches
    1: affine 4-space
    2: affine 4-space
    3: affine 4-space
  in the coordinate(s)
    1: [(x1//x0), (x2//x0), u, v]
    2: [(x0//x1), (x2//x1), u, v]
    3: [(x0//x2), (x1//x2), u, v]

julia> U = first(affine_charts(X))
Spectrum
  of multivariate polynomial ring in 4 variables (x1//x0), (x2//x0), u, v
    over rational field

julia> phi = homogenization_map(P, U);

julia> R = OO(U);

julia> phi.(gens(R))
4-element Vector{Tuple{MPolyDecRingElem{QQMPolyRingElem, AbstractAlgebra.Generic.MPoly{QQMPolyRingElem}}, MPolyDecRingElem{QQMPolyRingElem, AbstractAlgebra.Generic.MPoly{QQMPolyRingElem}}}}:
 (x1, x0)
 (x2, x0)
 (u, 1)
 (v, 1)
```
"""
function homogenization_map(P::AbsProjectiveScheme, U::AbsAffineScheme)
  # Projective schemes over a Field or ZZ or similar
  cache = _homogenization_cache(P)
  if haskey(cache, U)
    return cache[U]
  end
  error("patch not found or homogenization map not set")
end

function _homogenization_map(P::AbsProjectiveScheme, U::AbsAffineScheme, i::Int)
  # Determine those variables which come from the homogeneous
  # coordinates
  S = homogeneous_coordinate_ring(P)
  n = ngens(S)
  R = ambient_coordinate_ring(U)
  v = copy(gens(S))
  # prepare a vector of elements on which to evaluate the lifts
  popat!(v, i)  # remove the i-th variable
  function my_dehom(a::RingElem)
    parent(a) === OO(U) || error("element does not belong to the correct ring")
    p = lifted_numerator(a)
    q = lifted_denominator(a)
    deg_p = total_degree(p)
    deg_q = total_degree(q)
    deg_a = deg_p - deg_q
    ss = S[i] # the homogenization variable

    # preliminary lifts, not yet homogenized!
    pp = lift(evaluate(p, v))
    qq = lift(evaluate(q, v))
    # homogenize numerator and denominator
    pp = sum([c*m*ss^(deg_p - total_degree(m)) for (c, m) in zip(coefficients(pp), monomials(pp))])
    qq = sum([c*m*ss^(deg_q - total_degree(m)) for (c, m) in zip(coefficients(qq), monomials(qq))])

    if deg_a > 0
      return (pp, qq*ss^deg_a)
    elseif deg_a <0
      return (pp * ss^(-deg_a), qq)
    end
    return (pp, qq)
  end
  return my_dehom
end

function _homogenization_map(P::AbsProjectiveScheme{<:MPolyAnyRing, <:MPolyDecRing}, U::AbsAffineScheme, i::Int)
  # Determine those variables which come from the homogeneous
  # coordinates
  S = homogeneous_coordinate_ring(P)
  n = ngens(S)
  R = ambient_coordinate_ring(U)
  x = gens(R)
  s = x[1:n-1]
  x = x[n:end]
  B = base_ring(P)
  y = gens(B)
  t = gens(S)

  w = vcat([1 for j in 1:n-1], [0 for j in n:ngens(R)])
  v = copy(gens(S))
  # prepare a vector of elements on which to evaluate the lifts
  popat!(v, i)
  v = vcat(v, S.(gens(B)))
  function my_dehom(a::RingElem)
    parent(a) === OO(U) || error("element does not belong to the correct ring")
    p = lifted_numerator(a)
    q = lifted_denominator(a)
    deg_p = weighted_degree(p, w)
    deg_q = weighted_degree(q, w)
    deg_a = deg_p - deg_q
    ss = S[i] # the homogenization variable

    # preliminary lifts, not yet homogenized!
    pp = evaluate(p, v)
    qq = evaluate(q, v)

    # homogenize numerator and denominator
    pp = sum([c*m*ss^(deg_p - total_degree(m)) for (c, m) in zip(coefficients(pp), monomials(pp))])
    qq = sum([c*m*ss^(deg_q - total_degree(m)) for (c, m) in zip(coefficients(qq), monomials(qq))])

    if deg_a > 0
      return (pp, qq*ss^deg_a)
    elseif deg_a <0
      return (pp * ss^(-deg_a), qq)
    end
    return (pp, qq)
  end
  return my_dehom
end

function _homogenization_map(P::AbsProjectiveScheme{<:MPolyAnyRing, <:MPolyQuoRing}, U::AbsAffineScheme, i::Int)
  # Determine those variables which come from the homogeneous
  # coordinates
  S = homogeneous_coordinate_ring(P)
  n = ngens(S)
  R = ambient_coordinate_ring(U)
  x = gens(R)
  s = x[1:n-1]
  x = x[n:end]
  B = base_ring(P)
  y = gens(B)
  t = gens(S)

  w = vcat([1 for j in 1:n-1], [0 for j in n:ngens(R)])
  v = copy(gens(S))
  # prepare a vector of elements on which to evaluate the lifts
  popat!(v, i)
  v = vcat(v, S.(gens(B)))
  function my_dehom(a::RingElem)
    parent(a) === OO(U) || error("element does not belong to the correct ring")
    p = lifted_numerator(a)
    q = lifted_denominator(a)
    deg_p = weighted_degree(p, w)
    deg_q = weighted_degree(q, w)
    deg_a = deg_p - deg_q
    ss = S[i] # the homogenization variable

    # preliminary lifts, not yet homogenized!
    pp = evaluate(p, v)
    qq = evaluate(q, v)

    # homogenize numerator and denominator
    pp = sum([c*m*ss^(deg_p - total_degree(m)) for (c, m) in zip(coefficients(lift(pp)), monomials(lift(pp)))])
    qq = sum([c*m*ss^(deg_q - total_degree(m)) for (c, m) in zip(coefficients(lift(qq)), monomials(lift(qq)))])

    if deg_a > 0
      return (pp, qq*ss^deg_a)
    elseif deg_a <0
      return (pp * ss^(-deg_a), qq)
    end
    return (pp, qq)
  end
  return my_dehom
end

function getindex(X::AbsProjectiveScheme, U::AbsAffineScheme)
  Xcov = covered_scheme(X)
  for C in coverings(Xcov)
    for j in 1:n_patches(C)
      if U === C[j]
        return C, j
      end
    end
  end
  return nothing, 0
end

# comparison of projective spaces
function ==(X::AbsProjectiveScheme{<:Any,<:MPolyDecRing}, Y::AbsProjectiveScheme{<:Any,<:MPolyDecRing})
  return homogeneous_coordinate_ring(X) === homogeneous_coordinate_ring(Y)
end

# comparison of subschemes of projective space
function ==(X::AbsProjectiveScheme, Y::AbsProjectiveScheme)
  ambient_space(X) == ambient_space(Y) || return false
  IX = defining_ideal(X)
  IY = defining_ideal(Y)
  R = homogeneous_coordinate_ring(ambient_space(X))
  irrelevant_ideal = ideal(R,gens(R))
  IXsat = saturation(IX, irrelevant_ideal)
  IYsat = saturation(IY, irrelevant_ideal)
  return IXsat == IYsat
end

function is_subscheme(X::AbsProjectiveScheme, Y::AbsProjectiveScheme)
  ambient_space(X) == ambient_space(Y) || return false
  IX = defining_ideal(X)
  IY = defining_ideal(Y)
  R = homogeneous_coordinate_ring(ambient_space(X))
  irrelevant_ideal = ideal(R,gens(R))
  IXsat = saturation(IX, irrelevant_ideal)
  IYsat = saturation(IX, irrelevant_ideal)
  return issubset(IYsat, IXsat)
end

function Base.intersect(X::AbsProjectiveScheme, Y::AbsProjectiveScheme)
  return intersect([X, Y])
end

function Base.intersect(comp::Vector{<:AbsProjectiveScheme})
  @assert length(comp) > 0 "list of schemes must not be empty"
  IP = ambient_space(first(comp))
  @assert all(x->ambient_space(x)==IP, comp[2:end]) "schemes must have the same ambient space"
  S = homogeneous_coordinate_ring(IP)
  I = sum([defining_ideal(x) for x in comp])
  result = subscheme(IP, I)
  set_attribute!(result, :ambient_space, IP)
  return result
end

function Base.union(X::AbsProjectiveScheme, Y::AbsProjectiveScheme)
  return union([X, Y])
end

function Base.union(comp::Vector{<:AbsProjectiveScheme})
  @assert length(comp) > 0 "list of schemes must not be empty"
  IP = ambient_space(first(comp))
  @assert all(x->ambient_space(x)==IP, comp[2:end]) "schemes must have the same ambient space"
  S = homogeneous_coordinate_ring(IP)
  I = intersect([defining_ideal(x) for x in comp])
  result = subscheme(IP, I)
  set_attribute!(result, :ambient_space, IP)
  return result
end

function Base.hash(X::AbsProjectiveScheme, u::UInt)
  return u
end


