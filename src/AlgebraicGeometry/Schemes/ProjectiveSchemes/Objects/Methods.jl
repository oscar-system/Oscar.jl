################################################################################
# Printing and IO
################################################################################

function Base.show(io::IO, ::MIME"text/plain", P::AbsProjectiveScheme{<:Any, <:MPolyQuoRing})
  println(io, "Projective scheme")  # at least one new line is needed
  println(io, "  over ", base_ring(P))
  print(io, "  defined by ")
  print(io, defining_ideal(P)) # the last print statement must not add a new line
end

function Base.show(io::IO, P::AbsProjectiveScheme{<:Any, <:MPolyQuoRing})
  if get(io, :supercompact, false)
    # no nested printing
    print(io, "Projective scheme")
  else
    # nested printing allowed, preferably supercompact
    print(io, "Projective scheme in ")
    print(IOContext(io, :supercompact => true), ambient_space(P), " over ", base_ring(P))
  end
end

# Projective space
function Base.show(io::IO, ::MIME"text/plain", P::AbsProjectiveScheme{<:Any, <:MPolyDecRing})
  println(io, "Projective space of dimension $(relative_ambient_dimension(P))")  # at least one new line is needed
  print(io, "  with homogeneous coordinates ")
  for x in homogeneous_coordinates(P)
    print(io, x, " ")
  end
  println(io, "")
  print(io, "  over ")
  print(io, base_ring(P)) # the last print statement must not add a new line
end

function Base.show(io::IO, P::AbsProjectiveScheme{<:Any, <:MPolyDecRing})
  if get(io, :supercompact, false) # no nested printing
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, "ℙ$(ltx["\\^$(relative_ambient_dimension(P))"])")
    else
      print(io, "IP^$(relative_ambient_dimension(P))")
    end
  else
    if is_unicode_allowed()
      ltx = Base.REPL_MODULE_REF.x.REPLCompletions.latex_symbols
      print(io, "ℙ")
      n = relative_ambient_dimension(P)
      for d in reverse(digits(n))
        print(io, ltx["\\^$d"])
      end
      print(io, " over ")
      print(IOContext(io, :supercompact => true), base_ring(P))
    else
      # nested printing allowed, preferably supercompact
      print(io, "Projective $(relative_ambient_dimension(P))-space over ")
      print(IOContext(io, :supercompact => true), base_ring(P))
    end
  end
end


@doc raw"""
    dehomogenization_map(X::AbsProjectiveScheme, U::AbsSpec)

Return the restriction morphism from the graded coordinate ring of ``X`` to `𝒪(U)`.

# Examples
```jldoctest
julia> P = projective_space(QQ, ["x0", "x1", "x2"])
Projective space of dimension 2
  with homogeneous coordinates x0 x1 x2
  over Rational field

julia> X = covered_scheme(P);

julia> U = first(affine_charts(X))
Spec of Quotient of multivariate polynomial ring by ideal with 0 generators

julia> phi = dehomogenization_map(P, U);

julia> S = homogeneous_coordinate_ring(P);

julia> phi(S[2])
(x1//x0)

```
"""
function dehomogenization_map(X::AbsProjectiveScheme, U::AbsSpec)
  cache = _dehomogenization_cache(X)
  if haskey(cache, U)
    return cache[U]
  end
  charts = affine_charts(covered_scheme(X))
  any(x->(x===U), charts) || error("second argument is not an affine chart of the first")
  i = findfirst(k->(charts[k] === U), 1:relative_ambient_dimension(X)+1) - 1
  S = homogeneous_coordinate_ring(X)
  C = default_covering(covered_scheme(X))
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  phi = hom(S, OO(U), s, check=false)
  cache[U] = phi
  return phi
end

function dehomogenization_map(
    X::AbsProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  i in 0:relative_ambient_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  cache = _dehomogenization_cache(X)
  if haskey(cache, U)
    return cache[U]
  end
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  phi = hom(S, OO(U), pullback(p[U]), s, check=false)
  cache[U] = phi
  return phi
end


function dehomogenization_map(
    X::AbsProjectiveScheme{CRT}, 
    U::AbsSpec
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  return dehomogenization_map(X, X[U][2]-1)
end

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


@doc raw"""
    homogenization_map(P::AbsProjectiveScheme, U::AbsSpec)

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
not an instance of `Hecke.Map`.

# Examples
```jldoctest
julia> A, _ = QQ["u", "v"];

julia> P = projective_space(A, ["x0", "x1", "x2"])
Projective space of dimension 2
  with homogeneous coordinates x0 x1 x2
  over Multivariate polynomial ring in 2 variables over QQ

julia> X = covered_scheme(P);


julia> U = first(affine_charts(X))
Spec of Localization of quotient of multivariate polynomial ring at products of 1 element

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
function homogenization_map(P::AbsProjectiveScheme, U::AbsSpec)
  error("method not implemented for this type of input")
end

function homogenization_map(P::AbsProjectiveScheme{<:Field}, U::AbsSpec)
  error("method not implemented for projective schemes over fields")
end

function homogenization_map(P::AbsProjectiveScheme{<:Any, <:MPolyDecRing}, U::AbsSpec)
  cache = _homogenization_cache(P)
  if haskey(cache, U)
    return cache[U]
  end
  # Find the chart where U belongs to
  X = covered_scheme(P)
  i = findfirst(V->(U===V), affine_charts(X))
  i === nothing && error("the given affine scheme is not one of the standard affine charts")

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
  v = gens(S)
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
  cache[U] = my_dehom
  return my_dehom
end

function homogenization_map(P::AbsProjectiveScheme{<:Any, <:MPolyQuoRing}, U::AbsSpec)
  cache = _homogenization_cache(P)
  if haskey(cache, U)
    return cache[U]
  end
  # Find the chart where U belongs to
  X = covered_scheme(P)
  i = findfirst(V->(U===V), affine_charts(X))
  i === nothing && error("the given affine scheme is not one of the standard affine charts")
  
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
  v = gens(S)
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
  cache[U] = my_dehom
  return my_dehom
end

function getindex(X::AbsProjectiveScheme, U::AbsSpec)
  Xcov = covered_scheme(X)
  for C in coverings(Xcov)
    for j in 1:npatches(C)
      if U === C[j] 
        return C, j
      end
    end
  end
  return nothing, 0
end


