@doc raw"""
    dehomogenization_map(X::AbsProjectiveScheme, U::AbsSpec)

Return the restriction morphism from the graded coordinate ring of ``X`` to `𝒪(U)`.
"""
function dehomogenization_map(
    X::AbsProjectiveScheme{<:Ring}, 
    U::AbsSpec
  )
  if !isdefined(X, :dehomogenization_cache)
    X.dehomogenization_cache = IdDict()
  end
  cache = X.dehomogenization_cache
  if haskey(cache, U)
    return cache[U]
  end
  charts = affine_charts(covered_scheme(X))
  any(x->(x===U), charts) || error("second argument is not an affine chart of the first")
  i = findfirst(k->(charts[k] === U), 1:relative_ambient_dimension(X)+1) - 1
  S = homogeneous_coordinate_ring(X)
  C = default_covering(covered_scheme(X))
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  phi = hom(S, OO(U), s)
  cache[U] = phi
  return phi
end

@doc raw"""
    dehomogenization_map(X::AbsProjectiveScheme, i::AbsSpec)

Return the restriction morphism from the graded coordinate ring of ``X`` to `𝒪(Uᵢ)`.
Where `Uᵢ` is the `i`-th affine chart of `X`.
"""
function dehomogenization_map(
    X::AbsProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  if !isdefined(X, :dehomogenization_cache)
    X.dehomogenization_cache = IdDict()
  end
  i in 0:relative_ambient_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  cache = X.dehomogenization_cache
  if haskey(cache, U)
    return cache[U]
  end
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  phi = hom(S, OO(U), pullback(p[U]), s)
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

function dehomogenization_map(
    X::AbsProjectiveScheme{CRT},
    i::Int
  ) where {
    CRT<:AbstractAlgebra.Ring
  }
  i in 0:relative_ambient_dimension(X) || error("the given integer is not in the admissible range")
  S = homogeneous_coordinate_ring(X)
  C = default_covering(covered_scheme(X))
  U = C[i+1]
  if !isdefined(X, :dehomogenization_cache)
    X.dehomogenization_cache = IdDict()
  end
  cache = X.dehomogenization_cache
  if haskey(cache, U)
    return cache[U]
  end
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  phi = hom(S, OO(U), s)
  cache[U] = phi
  return phi
end


@doc raw"""
    homogenization_map(P::AbsProjectiveScheme, U::AbsSpec) -> function

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
"""
function homogenization_map(P::AbsProjectiveScheme{<:Any, <:MPolyDecRing}, U::AbsSpec)
  if !isdefined(P, :homogenization_cache)
    P.homogenization_cache = IdDict()
  end
  cache = P.homogenization_cache
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
  if !isdefined(P, :homogenization_cache)
    P.homogenization_cache = IdDict()
  end
  cache = P.homogenization_cache
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

function getindex(X::ProjectiveScheme, U::AbsSpec)
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


