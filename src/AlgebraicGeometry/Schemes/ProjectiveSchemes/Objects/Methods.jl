function dehomogenize(
    X::AbsProjectiveScheme{<:Ring}, 
    U::AbsSpec
  )
  charts = affine_charts(covered_scheme(X))
  any(x->(x===U), charts) || error("second argument is not an affine chart of the first")
  i = findfirst(k->(charts[k] === U), 1:relative_ambient_dimension(X)+1) - 1
  S = graded_coordinate_ring(X)
  C = default_covering(covered_scheme(X))
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  return hom(S, OO(U), s)
end

function dehomogenize(
    X::AbsProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  i in 0:relative_ambient_dimension(X) || error("the given integer is not in the admissible range")
  S = graded_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  return hom(S, OO(U), pullback(p[U]), s)
end

function dehomogenize(
    X::AbsProjectiveScheme{CRT}, 
    U::AbsSpec
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  return dehomogenize(X, X[U][2]-1)
end

function dehomogenize(
    X::AbsProjectiveScheme{CRT},
    i::Int
  ) where {
    CRT<:AbstractAlgebra.Ring
  }
  i in 0:relative_ambient_dimension(X) || error("the given integer is not in the admissible range")
  S = graded_coordinate_ring(X)
  C = default_covering(covered_scheme(X))
  U = C[i+1]
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  return hom(S, OO(U), s)
end


@doc Markdown.doc"""
    homogenize(P::AbsProjectiveScheme, U::AbsSpec)

Given an affine chart ``U ⊂ P`` of an `AbsProjectiveScheme` 
``P``, return a method ``h`` for the homogenization of elements 
``a ∈ 𝒪(U)``. 

That means ``h(a)`` returns a pair ``(p, q)`` representing a fraction 
``p/q ∈ S`` of the `ambient_coordinate_ring` of ``P`` such 
that ``a`` is the dehomogenization of ``p/q``.

**Note:** For the time being, this only works for affine 
charts which are of the standard form ``sᵢ ≠ 0`` for ``sᵢ∈ S``
one of the homogeneous coordinates of ``P``.

**Note:** Since this map returns representatives only, it 
is not a mathematical morphism and, hence, in particular 
not an instance of `Hecke.Map`.

**Note:** Since `fraction_field` relies on some implementation 
of division for the elements, we can not return the fraction 
directly. 
"""
function homogenize(P::AbsProjectiveScheme{<:Any, <:MPolyDecRing}, U::AbsSpec)
  # TODO: Ideally, one needs to provide this function 
  # only once for every pair (P, U), so we should think of 
  # some internal way for caching. The @attr macro is not 
  # suitable for this, because it is not sensitive for 
  # which U is put in. 

  # Find the chart where a belongs to
  X = covered_scheme(P)
  i = findfirst(V->(U===V), affine_charts(X))
  i === nothing && error("the given affine scheme is not one of the standard affine charts")
  
  # Determine those variables which come from the homogeneous 
  # coordinates
  S = graded_coordinate_ring(P)
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
  return my_dehom
end

function homogenize(P::AbsProjectiveScheme{<:Any, <:MPolyQuoRing}, U::AbsSpec)
  # TODO: Ideally, one needs to provide this function 
  # only once for every pair (P, U), so we should think of 
  # some internal way for caching. The @attr macro is not 
  # suitable for this, because it is not sensitive for 
  # which U is put in. 

  # Find the chart where a belongs to
  X = covered_scheme(P)
  i = findfirst(V->(U===V), affine_charts(X))
  i === nothing && error("the given affine scheme is not one of the standard affine charts")
  
  # Determine those variables which come from the homogeneous 
  # coordinates
  S = graded_coordinate_ring(P)
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

#function dehomogenize(
#    X::ProjectiveScheme{CRT}, 
#    U::AbsSpec
#  ) where {
#    CRT<:MPolyQuoLocRing
#  }
#  # look up U in the coverings of X
#  cover_of_U, index_of_U = X[U]
#  Xcov = covered_scheme(X)
#  S = graded_coordinate_ring(X)
# 
#  s = Vector{elem_type(OO(U))}()
#  if cover_of_U === standard_covering(X)
#    S = ambient_coordinate_ring(X)
#    C = standard_covering(X)
#    p = covered_projection_to_base(X)
#    s = vcat(gens(OO(U))[1:index_of_U-1], [one(OO(U))], gens(OO(U))[index_of_U:relative_ambient_dimension(X)])
#    return hom(S, OO(U), pullback(p[U]), s)
#  else
#    ref = Xcov[cover_of_U, standard_covering(X)]
#    V = codomain(ref[U])
#    index_of_V = standard_covering(X)[V]
#    t = vcat(gens(OO(V))[1:index_of_V-1], [one(OO(V))], gens(OO(V))[index_of_V:relative_ambient_dimension(X)])
#    s = pullback(ref[U]).(t)
#    pb = compose(ref[U], covered_projection_to_base(X)[V])
#    return hom(S, OO(U), pullback(pb), s)
#  end
#end

